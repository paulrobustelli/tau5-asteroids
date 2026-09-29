import itertools
import numpy as np
from selection import Objective,select
rng=np.random.default_rng(21)
x=rng.normal(size=(12,5));curves=rng.uniform(1,3,(12,7));target=np.array([1,4,8]);y=2.3*curves[target].mean(axis=0)
obj=Objective(x,x[target].mean(axis=0),np.ones(5),curves,y,np.ones(7))
assert abs(obj.evaluate(obj.features[target].mean(axis=0))[0])<1e-25
all_sets=np.array(list(itertools.combinations(range(12),3)))
vals=obj.evaluate(np.array([obj.features[z].mean(axis=0) for z in all_sets]))
best,trace=select(obj,3,seed=8,population=30,generations=100)
assert len(set(best))==3
assert np.all(np.diff(trace)<=1e-10)
assert np.isclose(trace[-1],vals.min(),atol=1e-12),(best,trace[-1],vals.min())
a,_=select(obj,12,seed=1,population=4,generations=2);assert len(a)==12
try:select(obj,13)
except ValueError:pass
else:raise AssertionError('oversized ensemble accepted')
print('PASS: known-target exhaustive optimum, uniform SAXS scaling, unique cardinality, monotonic elite, and full-pool edge case')
