# Ideas about how to optimize the process
First, before spending a lot of time on optimizing something, I shall profile and identify expensive computations.
## Data loading / format
There's probably a lot to improve in the way I load data. Creating this array of dicts of dicts seems very inefficient. Data like constants are copied around. Some computations could be already pre-carried at this point, like adding all the global variables. This would reduce the space and time cost of moving all of this data around.
Also, the way we do batches seems very inefficient. Perhaps I should use a proper data loader library, doing that in a more sensible way. 
## Pipeline
We could realise that a ctx is not relevant before/without making a pipeline forward pass for this ctx. Akash gave me some formula, I could use it, or do something else.