# Config streamlining
in the rl_project, we used:
* hydra - not sure what it does yet
* omegaconf (from omegaconf import DictConfig, OmegaConf)
* wandb for logging and plotting
proposed config structure:
model:
    n_hidden
    hidden_size # Perhaps replace with an array describing layers?
    batchNorm
    activation?other?
data: (does it need to be nested?)
    site(s)
    base_path
    nPoints
pipeline:
    modelMapping? default=softplus
    which pipeline?
    use Vmax?
train:
    lr
    epochs
    weight_decay
    batch_size
    loss?opt?iterator?