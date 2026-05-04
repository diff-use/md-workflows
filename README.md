# md-workflows

This project includes a `Dockerfile` and a Python CLI entry point for running the full MD workflow pipeline.

## 1) Build the Docker image

From the project root (where the `Dockerfile` is):

```bash
docker pull diffuseproject/md:0.0.1 .
```

## 2) Start a container

Run the container interactively, mounting the current project directory so inputs/outputs are available on your host:

```bash
docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  --gpus all \
  --name md_container \
  -e HOME=/workspace \
  -v "$(pwd):/workspace" \
  -w /workspace \
  diffuseproject/md:0.0.1 \
  bash
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
md_workflows.mdmx \
  --param-pdb-id 6B8X \
  --ix 1 \
  --ntomp 26 \
  --resolv-ntmpi 8 \
  --resolv-ntomp 1
```

To see all available flags:

```bash
md_workflows.mdmx --help
```
