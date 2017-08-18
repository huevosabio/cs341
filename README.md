# CS 341

# TODO:

- *Controller*: Currently, demand that is _not_ satisfied simply gets ignored. Moreover, we ignore `FlowsOut`, but not `FlowsIn`. When we "ignore" a demand, we should be postponing it, such that the corresponding demand get's postponed. This might ruin our integrality either by forcing us to use non-integral predictions or making our LP non-TUM. Also, it would make the problem significantly larger since we will likely end up with at least `M+1` classes, where `M` is the number of stations.

- *Controller*: Many travel requests are within the same district. It was necessary to scrap the multi-commodity flow approach 


- *Predictor*: Currently, the predictor is being trained _only_ on demand by (time,station). However, the controller depends on (time, origin, destination, travel_time). Alternatively, the option might be to simply predict the FlowsIn and FlowsOut independently. This would, however, not pair at all with the idea of being able to "postpone" demand. 

## Strategies



|Controller \ Predictor   | TO | TOD | TODTT | FlowsIn/Out | 
| ----------------------- |:--:|:---:|:-----:|:-----------:|
|**Origin Slack**         | E  |     |       |             |
|**TOD Postpone**         |    |     |       |             |
|**Ignore Flows-In**      | E  |     |       |             |



## Notes

_08/16/2017_

The multi-commodity approach for the controller was not possible, since passengers could simply "wait" and that would be considered as if the demand was satisfied. Thus, the approach to model waiting is now to use a specific set of variables that keeps track of waiting passengers. Not surprisingly, this *i) makes the problem slower* and *ii) breaks total unimodularity*. Changing from waiting to droping the customers solves both i) and ii). This, however, breaks our hope to be able to keep track of waiting customers. I will attempt now to simply provide waiting _only_ for the passengers that are already there. 

The results for the approach with waiting only for the immediate customers appear promising: it does not suffer from either i) or ii), and in general avoids waiting as much as possible. The number of dropped customers is, however, worrying. For example, we are seeing serivce percentages of ~70% Poisson demands of lambda=0.2 (but ~95% with lambda=0.1). I will now try to implement this with the actual (o,d,t) demands (travel time still uncertain). 