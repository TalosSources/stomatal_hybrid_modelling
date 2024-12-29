* ask about the missing data thing. timesteps, contexts. how to determine which contexts are valid? 
* ask about the daily averaging. how to do it? average the predictors as well? seems quite sketchy.
    a way to do it could be to, for all the timesteps of a full day (or other period), compute a Q_LE,
    then average all those output Q_LE to obtain the model's estimation for the day's average,
    and then backward on that, using the error w.r.t the average of y Q_LE for that day.