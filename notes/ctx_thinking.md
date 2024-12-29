I don't yet fully understand how things work.
But let's make the most general assumption:
Each timestep may contain valid or invalid data for each context independantly.
Ideally, it would be best to have some ground truth flag indicating which ctx is valid along with the actual data.
With this, I could safely discard other contexts, and ignore any value that makes the pipeline fail as an ill-value.
However, I don't think I have this.
TODO: I should ask Akash if I can determine this directly from the data.
But assuming I can't, we assume each time-step can have data for low and high vegetation, and that potentially, there could be both valid shaded, sunny component to the flux, or none.
The only way to figure what ctx is valid, is to look for 0. or NaN values (or potentially inf), in the timestep values.

# Mapping
Here are the variables depending on both sun/shd and H/L context:
* Cc
* IPAR
* Psi_L
* Vmax
* (rs)

Here are the variables depending only on H/L context:
* ra
* rb
* constants: Psi_sto_50, Psi_sto_100, CT, Ha, FI, Do, a1, go, gmes, rjv, DS

Some values can be problematic: 
* if any value is NaN, the corresponding context must be invalid, since all values are used. If all ctxs contain a NaN value, or a global value is Nan, perhaps in some cases it can be treated as 0, but in principle the whole timestep is invalid
* ra is 0, there's a division by 0 in the PM equation, so the ctx must be treated as invalid.
* if other values are 0, it's tricky to choose what to do. Perhaps some values are expected to be 0. But sometimes, it means the value is missing. In this case, using it as is can disturb the learning process.
* IPAR is vital to the pb module, can't be NaN -> implies the context is invalid

# Method
Let's decide what to do in general. 
Devise a list of checks for ctx dependant variables, s.t. if at least one fails, the context is considered invalid.
Devise a list of checks for a timestep (one of them is that at least 1 context is valid for this timestep), s.t. if at least one fails, the timestep is considered invalid.
If a timestep is valid, include in its context dict all the valid contexts, and add the timestep to the dataset.
That way, the rest of the pipeline, which is just treating ctx dicts in an agnostic way, shouldn't require any change.

# Big problem
I think it prevents us from doing batches the way we were doing it?
Suppose in batch $b$, there's a timestep with only ctx sun_h and another with only shd_L.
Then, the tensor of size |b| representing, for example, Cc_sun_H, can't have a meaningful entry
for the shd_L timestep. we could always do ugly tricks like put in a dummy nan entry, indicating this
timestep is invalid for this context, and then replace nan values in the end. However, this is very cumbersome,
ugly, and I don't even know how to deal with gradients in this case.
So let's keep things simple. We can do batches of size 1 at a time in the pipeline, but we loose the vectorized operation speed-up.
So instead, we can do this: go through several timesteps, agreggate them in some way?
but that's quite complicated, it requires to know which context-timestep came from where to continue the computation of Q_LE afterwards?
the knot of the problem is that we can't at the same time not compute the pipeline for a context and do it, if we have several contradicting timesteps. There's this branching that causes issues.
Here's an idea: we group the timesteps by the ctx set they have. There's at most 15 different ctx sets, and in practice, surely less than that.
Then, at batch sampling, we choose a set randomly, weighted by the set size. Then, we select a batch from this set. 
We can safely perform our vectorized operations that way.

# Questions that subsist
I've been weighing low and tall vegetation in the same way. Is it really how it should be done?


What to do if there's no data for a timestep? we get a gradient error. Ideally, we should discard it from the dataset.