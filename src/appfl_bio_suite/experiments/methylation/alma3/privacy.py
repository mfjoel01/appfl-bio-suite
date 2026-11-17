"""Full-batch Gaussian mechanism; fixed-size replace-one patient adjacency.

No subsampling amplification. Per-step sensitivity of the clipped gradient sum
is 2C. Gaussian noise std is 2C*sigma. rho = steps/(2*sigma**2).
Calibration histograms use disjoint correct/incorrect bins: sensitivity sqrt(2).
"""

import math
import secrets

import torch


def private_generator(device="cpu"):
    # Independent of published study seeds. This is a research PRNG, not a claim
    # of a hardened cryptographic Gaussian sampler.
    return torch.Generator(device=device).manual_seed(secrets.randbits(63))


def epsilon(rho, delta):
    if not math.isfinite(rho) or rho < 0 or not 0 < delta < 1:
        raise ValueError("invalid privacy accounting inputs")
    return rho + 2 * math.sqrt(rho * math.log(1 / delta))


def accounting(sigma, steps, calibration_queries, delta):
    if sigma == 0:
        return {"epsilon": None, "rho": None, "delta": delta, "private": False}
    # Training and calibration cohorts are disjoint: parallel composition.
    train_rho = steps / (2 * sigma**2)
    calibration_rho = calibration_queries / (2 * sigma**2)
    rho = max(train_rho, calibration_rho)
    return {
        "epsilon": epsilon(rho, delta),
        "rho": rho,
        "delta": delta,
        "training_rho": train_rho,
        "calibration_rho": calibration_rho,
        "private": True,
        "adjacency": "replace-one patient, fixed public site sizes",
        "scope": "model updates and calibration only; synthetic evaluation is public",
    }


def train(model, inputs, labels, *, steps, lr, clip, sigma, mu=0.0):
    """One full-data clipped-gradient step per iteration; optimizer reset per round.

    Clipping also applies to nonprivate arms, so noise is the only DP difference.
    The data-independent FedProx gradient is added after clipping/noising.
    """
    model.eval()  # no batch statistics or dropout interactions between patients
    parameters = [p for p in model.parameters() if p.requires_grad]
    anchor = [p.detach().clone() for p in parameters]
    optimizer = torch.optim.SGD(parameters, lr=lr)
    generator = private_generator(parameters[0].device)
    for _ in range(steps):
        sums = [torch.zeros_like(p) for p in parameters]
        for i in range(len(labels)):
            optimizer.zero_grad(set_to_none=True)
            logp = model(*(v[i : i + 1] for v in inputs))
            loss = -logp[0, labels[i]]
            if not torch.isfinite(loss):
                raise ValueError("nonfinite patient loss")
            loss.backward()
            gradients = [p.grad if p.grad is not None else torch.zeros_like(p) for p in parameters]
            norm = torch.stack([g.square().sum() for g in gradients]).sum().sqrt()
            factor = (clip / norm.clamp_min(1e-12)).clamp_max(1)
            for total, grad in zip(sums, gradients, strict=True):
                total.add_(grad * factor)
        optimizer.zero_grad(set_to_none=True)
        for p, total, initial in zip(parameters, sums, anchor, strict=True):
            if sigma:
                total.add_(
                    torch.randn(
                        total.shape, generator=generator, device=total.device, dtype=total.dtype
                    ),
                    alpha=2 * clip * sigma,
                )
            p.grad = total / len(labels) + mu * (p.detach() - initial)
        optimizer.step()
    if any(not torch.isfinite(p).all() for p in parameters):
        raise ValueError("nonfinite trained parameters")
