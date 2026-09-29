"""Fixed-cardinality, equal-weight genetic selection; independent reimplementation.

No continuous conformer weights and no predictor-error stopping threshold.
"""
import numpy as np

class Objective:
    def __init__(self, shifts, observed, shift_scales, saxs, intensity, errors, saxs_weight=1.):
        self.cs=(np.asarray(shifts,dtype=float)-observed)/shift_scales
        self.saxs=np.asarray(saxs,dtype=float)/errors
        self.y=np.asarray(intensity,dtype=float)/errors
        self.ncs=self.cs.shape[1]
        self.features=np.concatenate([self.cs,self.saxs],axis=1)
        self.saxs_weight=saxs_weight
        if not np.isfinite(self.features).all() or not np.isfinite(self.y).all():
            raise ValueError('Nonfinite observables')
    def evaluate(self,means):
        means=np.atleast_2d(means);z=means[:,:self.ncs];p=means[:,self.ncs:]
        denom=np.sum(p*p,axis=1)
        scale=np.maximum(0,(p@self.y)/np.maximum(denom,1e-300))
        return np.sum(z*z,axis=1)+self.saxs_weight*np.sum((scale[:,None]*p-self.y)**2,axis=1)


def select(objective,n_select,seed=1,population=100,generations=1000,callback=None):
    f=objective.features;n=len(f);rng=np.random.default_rng(seed)
    if not 1<=n_select<=n:raise ValueError('Insufficient candidates for distinct selection')
    if population<4:raise ValueError('Population must be >=4')
    pool=np.array([np.sort(rng.choice(n,n_select,replace=False)) for _ in range(population)])
    means=np.array([f[x].mean(axis=0) for x in pool]);scores=objective.evaluate(means)
    trace=[]
    def mutate(parent,mean,allowed):
        available=np.setdiff1d(allowed,parent,assume_unique=True)
        if not len(available):return parent.copy(),mean.copy()
        # Mix single replacements with larger moves, rather than a fixed mutation rate.
        count=min(len(available),n_select,int(rng.choice([1,2,5,max(1,n_select//20)])))
        positions=rng.choice(n_select,count,replace=False);added=rng.choice(available,count,replace=False)
        child=parent.copy();newmean=mean+(f[added].sum(axis=0)-f[parent[positions]].sum(axis=0))/n_select
        child[positions]=added
        return np.sort(child),newmean
    for gen in range(generations+1):
        best=int(np.argmin(scores));trace.append(float(scores[best]))
        if callback and (gen%25==0 or gen==generations):callback(gen,pool[best],scores[best],trace)
        if gen==generations:break
        children=[];childmeans=[];internal=np.unique(pool)
        for j in range(population):
            for allowed in (internal,np.arange(n)):
                child,mean=mutate(pool[j],means[j],allowed);children.append(child);childmeans.append(mean)
            other=int(rng.integers(population));union=np.union1d(pool[j],pool[other])
            child=np.sort(rng.choice(union,n_select,replace=False));children.append(child);childmeans.append(f[child].mean(axis=0))
            child=np.sort(rng.choice(n,n_select,replace=False));children.append(child);childmeans.append(f[child].mean(axis=0))
        children=np.array(children);childmeans=np.array(childmeans);childscores=objective.evaluate(childmeans)
        # Random disjoint tournaments among the 4P offspring, plus explicit elitism.
        groups=rng.permutation(len(children)).reshape(population,4)
        wins=groups[np.arange(population),np.argmin(childscores[groups],axis=1)]
        winner=int(np.argmin(childscores));wins[0]=winner
        nextpool=children[wins];nextmeans=childmeans[wins];nextscores=childscores[wins]
        if scores[best]<nextscores.max():
            k=int(np.argmax(nextscores));nextpool[k]=pool[best];nextmeans[k]=means[best];nextscores[k]=scores[best]
        pool,means,scores=nextpool,nextmeans,nextscores
        if gen%100==0:
            means=np.array([f[x].mean(axis=0) for x in pool]);scores=objective.evaluate(means)
    chosen=pool[np.argmin(scores)]
    assert len(np.unique(chosen))==n_select
    return chosen,trace
