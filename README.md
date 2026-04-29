# md-workflows

This project includes a `Dockerfile` and a Python CLI entry point for running the full MD workflow pipeline.

## 1) Build the Docker image

From the project root (where the `Dockerfile` is):

```bash
docker build -t md-workflows:latest .
```

## 2) Start a container

Run the container interactively, mounting the current project directory so inputs/outputs are available on your host:

```bash
docker run --rm -it \
  -v "$(pwd):/workspace" \
  -w /workspace \
  md-workflows:latest bash
```

## 3) Install this project inside the container

Inside the container shell:

```bash
python -m pip install -e .
```

This registers the CLI entry points from `pyproject.toml`, including `md_workflows.run_all`.

## 4) Run the full workflow command

Inside the container shell:

```bash
md_workflows.run_all
```

You can also pass options:

```bash
md_workflows.run_all \
  --param-pdb-id 6B8X \
  --ix 1 \
  --ntomp 26 \
  --resolv-ntmpi 8 \
  --resolv-ntomp 1
```

To see all available flags:

```bash
md_workflows.run_all --help
```
