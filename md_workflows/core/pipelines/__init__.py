"""Hardened multi-step pipelines composed from the individual steps.

Pipelines are plain function composition over ``core.steps`` — no orchestration
framework. The same composition can later be wrapped as a Prefect ``@flow`` unchanged.
"""
