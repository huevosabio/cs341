# CS 341

# TODO:

- *Predictor*: Currently, the predictor is being trained _only_ on demand by (time,station). However, the controller depends on (time, origin, destination, travel_time). Alternatively, the option might be to simply predict the FlowsIn and FlowsOut independently. This would, however, not pair at all with the idea of being able to "postpone" demand. 

## Strategies



|Controller \ Predictor   | TO | TOD | TODTT | FlowsIn/Out | 
| ----------------------- |:--:|:---:|:-----:|:-----------:|
|**Origin Slack**         | F  |     |       |             |
|**TOD Postpone**         | -  |     |       |             |
|**Ignore Flows-In**      | F  |     |       |             |
|**TOD Drop + Wait**      | S  |     |       |             |


## Notes

_08/16/2017_

The multi-commodity approach for the controller was not possible, since passengers could simply "wait" and that would be considered as if the demand was satisfied. Thus, the approach to model waiting is now to use a specific set of variables that keeps track of waiting passengers. Not surprisingly, this *i) makes the problem slower* and *ii) breaks total unimodularity*. Changing from waiting to droping the customers solves both i) and ii). This, however, breaks our hope to be able to keep track of waiting customers. I will attempt now to simply provide waiting _only_ for the passengers that are already there. 

The results for the approach with waiting only for the immediate customers appear promising: it does not suffer from either i) or ii), and in general avoids waiting as much as possible. The number of dropped customers is, however, worrying. For example, we are seeing serivce percentages of ~70% Poisson demands of lambda=0.2 (but ~95% with lambda=0.1). I will now try to implement this with the actual (o,d,t) demands (travel time still uncertain). 

_08/28/2017_

I took the average TOD trips for the weekdays and used that as a predictor to feed the controller. Sadly, the results show that it performs _better_ than our TommyLSTM. Upon some inspection via plotting, I suspect that the lack of quality on the later prediction steps has strong repercursions. It seems that by not being able to predict correctly the demand in remote stations with enough anticipation, the controller sends needed rebalancing vehicles a bit too late, incurring customer wait time. This is notable on both the spikes in rebalancing at ~5:30AM and ~13:00 for both "MPC Perfect" and "MPC Naive" but not "MPC LSTM". "MPC LSTM", instead, shows a similarly shaped spike a little later.


## Results

### TOD Predictions, Predicted Demand can be dropped, Waiting Customers must be satisfied

| Metric \ Predictor    | Perfect | TommyLSTM | SeasLSTM | Reactive | NoPred | Naive  |
| --------------------  |:------: |:---------:|:--------:|:--------:|:------:|:------:|
|**Mean wait time**     |   3.7   |   56.3    |   29.4   |  283.3   | 543.5  |  48.2  |
|**Median wait time**   |   0.0   |    0.0    |    0.0   |   72.0   | 414.0  |   0.0  |
|**Mean service time**  | 432.6   |  432.5    |  432.5   |  432.5   | 432.3  | 432.5  |
|**Median service time**| 330.0   |  330.0    |  330.0   |  330.0   | 330.0  | 330.0  |


## comments

Uncertainty Aware MPC Data-Driven Control of AMoD

(uncertainnty is the )

- 1: How can we accommodate surprises 
- 2: Including uncertainty
- 3: Experiments: Optimize Horizon (competing objectives long-horizon vs quality of predictions), trade short and long memory in the predictions, how do you shape the penalties in order to shave the peak in waiting times (not accounting for time windows, but we can penalize)

be careful with the prediction part 

## Tasks


- Build TOD LSTM predictor; Priority: LOW
-- Inputs/Outputs are now (T,O,D) tensors
-- Follow Yarin's dropout strategy closely
-- Attempt to use time/weekday as an input
-- Attempt to use Autoencoder to improve accuracy (see "Time Series Forecasting Based on Augmented Long Short-Term Memory")
-- Use Bayesian Optimization for parameter search
-- UPDATE: Actually, this is harder than I thought. TOD meand that we have 66^2 values per time step... making our LSTM _much_ larger.

- Emphasize in writing:
-- That the MILP does not grow with the number of vehicles
-- Intuition as to why does it run _so_ fast.

## Useful resources
- https://gallery.cortanaintelligence.com/Tutorial/Forecasting-Short-Time-Series-with-LSTM-Neural-Networks-2
- https://robjhyndman.com/hyndsight/transformations/
