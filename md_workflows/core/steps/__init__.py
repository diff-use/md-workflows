"""Individual MD prep steps.

Each module exposes the same contract:

* a ``*Inputs`` model (explicit, typed file paths + ``workdir``),
* ``resolve_inputs(workdir, cfg)`` — map the conventional filenames in a working
  directory to explicit input paths,
* ``check_inputs(inputs)`` — guard that raises ``MissingInputError`` before any tool
  runs,
* the step function itself, which takes ``*Inputs`` and returns a ``StepResult``.
"""
