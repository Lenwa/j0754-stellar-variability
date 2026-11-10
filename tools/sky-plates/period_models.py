"""Fixed-period, camera-offset and temporal-replay tools; no period FAP claim."""
from collections import defaultdict
import numpy as np

PERIOD = 912.3158778
EPOCH = 2458115.99394


def bin_series(rows, width=30):
    if len({r["band"] for r in rows}) > 1:
        raise ValueError("Different passbands require separate models")
    daily = defaultdict(list)
    for r in rows:
        daily[r["camera"],r["day"]].append(r)
    groups = defaultdict(list)
    for (camera,_), data in daily.items():
        t = float(np.median([r["jd"] for r in data]))
        m = float(np.median([r["mag"] for r in data]))
        groups[camera,int(np.floor((t-EPOCH)/width))].append((t,m,len(data)))
    result = [dict(camera=camera,bin_index=i,jd=float(np.median([x[0] for x in data])),
                   mag=float(np.median([x[1] for x in data])),days=len(data),frames=sum(x[2] for x in data))
              for (camera,i),data in groups.items()]
    return sorted(result,key=lambda r:(r["jd"],r["camera"]))


def design(t,camera,kind,period=PERIOD,harmonics=1,knots=None):
    t = np.asarray(t,float);camera = np.asarray(camera)
    columns = [(camera==c).astype(float) for c in sorted(set(camera))]
    if kind == "periodic":
        if not np.isfinite(period) or period <= 0:
            raise ValueError("Period must be finite and positive")
        phase = 2*np.pi*(t-EPOCH)/period
        for h in range(1,harmonics+1):
            columns += [np.sin(h*phase),np.cos(h*phase)]
    elif kind == "slow":
        knots = np.linspace(t.min(),t.max(),5) if knots is None else np.asarray(knots)
        for i in range(1,len(knots)):
            columns.append(np.interp(t,knots,np.eye(len(knots))[i]))
    elif kind != "null":
        raise ValueError("Unknown model")
    return np.column_stack(columns)


def fit_model(t,y,camera,kind,period=PERIOD,harmonics=1,train=None,knots=None):
    t=np.asarray(t,float);y=np.asarray(y,float);camera=np.asarray(camera)
    train=np.ones(len(t),bool) if train is None else np.asarray(train,bool)
    x=design(t,camera,kind,period,harmonics,knots)
    beta,_,rank,_=np.linalg.lstsq(x[train],y[train],rcond=None)
    if rank < x.shape[1]:
        raise ValueError("Training data do not identify all camera offsets/model terms")
    return dict(beta=beta,prediction=x@beta,rank=int(rank),parameters=x.shape[1])


def scores(y,prediction,baseline):
    y=np.asarray(y);error=y-np.asarray(prediction);null=y-np.asarray(baseline)
    mse=float(np.mean(error**2));null_mse=float(np.mean(null**2))
    return dict(rmse=float(np.sqrt(mse)),null_rmse=float(np.sqrt(null_mse)),
                r2=1-mse/null_mse if null_mse else None,
                rmse_gain=1-np.sqrt(mse/null_mse) if null_mse else None)


def forecast(t,y,camera,kind,cut,period=PERIOD,harmonics=1):
    train=np.asarray(t)<cut
    model=fit_model(t,y,camera,kind,period,harmonics,train)
    null=fit_model(t,y,camera,"null",train=train)
    return dict(scores(y[~train],model["prediction"][~train],null["prediction"][~train]),
                training_bins=int(train.sum()),test_bins=int((~train).sum()),
                beta=model["beta"].tolist())


def blocked_cv(t,y,camera,kind,period=PERIOD,harmonics=1):
    groups=np.floor((np.asarray(t)-EPOCH)/90).astype(int)
    prediction=np.full(len(t),np.nan);baseline=np.full(len(t),np.nan);skipped=[]
    for group in sorted(set(groups)):
        held=groups==group
        try:
            model=fit_model(t,y,camera,kind,period,harmonics,~held)
            null=fit_model(t,y,camera,"null",train=~held)
        except ValueError as e:
            skipped.append(dict(block=int(group),reason=str(e)));continue
        prediction[held]=model["prediction"][held];baseline[held]=null["prediction"][held]
    valid=np.isfinite(prediction)
    if not valid.any():
        raise ValueError("No identifiable blocked predictions")
    return dict(scores(y[valid],prediction[valid],baseline[valid]),
                bins_scored=int(valid.sum()),blocks=len(set(groups)),skipped_blocks=skipped)


def replay_pairs(t,y,camera,period=PERIOD,tolerance=15):
    t=np.asarray(t);y=np.asarray(y);camera=np.asarray(camera)
    candidates=[]
    for c in sorted(set(camera)):
        indexes=np.flatnonzero(camera==c)
        for i in indexes:
            for j in indexes[np.abs(t[indexes]-t[i]-2*period)<=tolerance]:
                if t[j]>t[i]:
                    candidates.append((abs(t[j]-t[i]-2*period),int(i),int(j)))
    used=set();pairs=[]
    for distance,i,j in sorted(candidates):
        if i in used or j in used:
            continue
        used.update([i,j])
        pairs.append(dict(camera=str(camera[i]),first_index=i,last_index=j,first_jd=float(t[i]),
                          last_jd=float(t[j]),phase_mismatch_days=float(distance),delta=float(y[j]-y[i])))
    return pairs
