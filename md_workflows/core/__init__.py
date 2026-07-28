"""Pure, orchestrator-agnostic MD workflow logic.

Everything under ``core`` is a plain, typed Python function with no CLI parsing and
no orchestration-framework dependency. The CLI, the SDK, and (later) a Prefect layer
are all just different callers of these functions.
"""
