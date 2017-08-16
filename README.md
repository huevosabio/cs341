# CS 341

# TODO:

- *Controller*: Currently, demand that is _not_ satisfied simply gets ignored. Moreover, we ignore `FlowsOut`, but not `FlowsIn`. When we "ignore" a demand, we should be postponing it, such that the corresponding demand get's postponed. This might ruin our integrality either by forcing us to use non-integral predictions or making our LP non-TUM. Also, it would make the problem significantly larger since we will likely end up with at least `M+1` classes, where `M` is the number of stations.


- *Predictor*: Currently, the predictor is being trained _only_ on demand by (time,station). However, the controller depends on (time, origin, destination, travel_time). Alternatively, the option might be to simply predict the FlowsIn and FlowsOut independently. This would, however, not pair at all with the idea of being able to "postpone" demand. 

## Strategies



|Controller \ Predictor   | TO | TOD | TODTT | FlowsIn/Out | 
| ----------------------- |:--:|:---:|:-----:|:-----------:|
|**Origin Slack**         | E  |     |       |             |
|**TOD Postpone**         |    |     |       |             |
|**Ignore Flows-In**      | E  |     |       |             |
