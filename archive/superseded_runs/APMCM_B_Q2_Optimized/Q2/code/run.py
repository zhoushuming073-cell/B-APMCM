from reconstruct import ROOT,TAB,pd,reconstruct
from delay import delay_study,local_profiles

if __name__=='__main__':
    TAB.mkdir(exist_ok=True)
    frame=reconstruct(pd.read_csv(ROOT/'data/process.csv'))
    delay_study(frame)
    local_profiles(frame)
    from global_lag import full_profiles
    full_profiles(frame)
    from checks import checks
    checks(frame)
    from figures import figures,dose_response
    figures()
    dose_response()
