import json,statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from study import ROOT

def main():
    rows=json.loads((ROOT/'results/cancellation.json').read_text())
    d={(r['name'],r['policy']):r for r in rows}
    summaries=[]
    for fam in ['delay','multidelay']:
        names=[r['name'] for r in rows if r['family']==fam and r['policy']=='jet']
        summary=dict(family=fam,n=len(names),
            fixed=statistics.median(d[n,'oracle9_2']['panels'] for n in names),
            real=statistics.median(d[n,'local2']['panels'] for n in names),
            jet=statistics.median(d[n,'jet']['panels'] for n in names),
            panel_ratio=statistics.median(d[n,'oracle9_2']['panels']/d[n,'jet']['panels'] for n in names),
            replay_ratio=statistics.median(d[n,'oracle9_2']['replay_seconds']/d[n,'jet']['replay_seconds'] for n in names),
            construction_ratio=statistics.median(d[n,'raw2']['seconds']/d[n,'jet']['seconds'] for n in names),
            wins=sum(d[n,'jet']['panels']<d[n,'oracle9_2']['panels'] for n in names),
            ties=sum(d[n,'jet']['panels']==d[n,'oracle9_2']['panels'] for n in names))
        summaries.append(summary)
    (ROOT/'results/cancellation_summary.json').write_text(json.dumps(summaries,indent=2))
    plt.rcParams.update({'font.family':'serif','font.size':10,'pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(6.5,3.2))
    for fam,color,marker in [('delay','#526aa1','o'),('multidelay','#057a75','s')]:
        rr=[r for r in rows if r['family']==fam and r['policy']=='jet']
        ax.scatter([d[r['name'],'oracle9_2']['panels'] for r in rr],
            [r['panels'] for r in rr],color=color,marker=marker,label=fam,s=30,alpha=.8)
    ax.plot([0,150],[0,150],'k--',lw=.8)
    ax.set(xlabel='Fixed-9 panels, second-order enclosure',ylabel='Jet panels, second-order enclosure')
    ax.spines[['top','right']].set_visible(False)
    ax.legend(frameon=False);fig.tight_layout()
    fig.savefig(ROOT/'figures/cancellation.pdf',bbox_inches='tight')
    fig.savefig(ROOT/'figures/cancellation.png',dpi=180,bbox_inches='tight')
    table=['\\begin{tabular}{lrrrrrr}\\toprule',
           'Family & $n$ & Fixed-9 & Real local & Jet & Panel ratio & Replay ratio\\\\\\midrule']
    for r in summaries:
        table.append(f"{r['family'].capitalize()} & {r['n']} & {r['fixed']:g} & {r['real']:g} & {r['jet']:g} & {r['panel_ratio']:.2f} & {r['replay_ratio']:.2f}\\\\")
    table+=['\\bottomrule\\end{tabular}']
    (ROOT/'paper/cancellation_table.tex').write_text('\n'.join(table))
    print(json.dumps(summaries,indent=2))

if __name__=='__main__':main()
