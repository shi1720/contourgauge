"""All manuscript figures and numerical tables derive from saved results."""
import json,statistics,math
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from study import ROOT

plt.rcParams.update({'font.family':'serif','font.size':10,'axes.spines.top':False,
 'axes.spines.right':False,'savefig.bbox':'tight','pdf.fonttype':42})
COLORS={'raw':'#8c8c8c','anchor':'#d88732','global':'#526aa1','local':'#057a75'}
LABEL={'raw':'No gauge','anchor':'Anchor','global':'Midrange','local':'Local'}

def save(fig,name):
    fig.savefig(ROOT/'figures'/f'{name}.pdf')
    fig.savefig(ROOT/'figures'/f'{name}.png',dpi=180)
    plt.close(fig)

def main():
    rows=json.loads((ROOT/'results/study.json').read_text())
    oracle=json.loads((ROOT/'results/global_sweep.json').read_text())
    replay=json.loads((ROOT/'results/replay_timing.json').read_text())
    by={(r['name'],r['policy']):r for r in rows}
    # 1: geometric explanation; deliberately uses a certified panel.
    fig,ax=plt.subplots(1,2,figsize=(9,3.3))
    t=np.linspace(-1,1,1000)
    f=np.exp(5j*t)*(1+0.15*t+0.12j*t*t)
    g=1+0.15*t+0.12j*t*t
    ax[0].plot(f.real,f.imag,color=COLORS['raw'],lw=1.5)
    ax[0].set_title('(a) Original panel image')
    ax[1].plot(g.real,g.imag,color=COLORS['local'],lw=2)
    ax[1].add_patch(plt.Circle((1,0),.27,fill=False,color=COLORS['local'],ls='--'))
    ax[1].set_title('(b) After removing $e^{5it}$')
    for a in ax:
        a.plot(0,0,'k+',ms=8);a.set_aspect('equal');a.set_xlabel('Real part');a.set_ylabel('Imaginary part')
    ax[1].set_xlim(-.15,1.45);ax[1].set_ylim(-.5,.5)
    fig.tight_layout();save(fig,'mechanism')
    # 2: common exponent shift, lower bound illustrates only ungauged disks.
    fig,ax=plt.subplots(figsize=(6.5,3.4))
    for p in COLORS:
        rr=[r for r in rows if r['family']=='shift' and r['policy']==p]
        ax.plot([r['shift']+1 for r in rr],[r['panels'] for r in rr],'o-',label=LABEL[p],color=COLORS[p])
    ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Common exponent shift $\\kappa+1$')
    ax.set_ylabel('Certified contour panels');ax.legend(ncol=2,frameon=False)
    ax.grid(alpha=.15);fig.tight_layout();save(fig,'shift')
    # 3: paired comparison with the offline best of nine fixed centers.
    fig,axes=plt.subplots(1,2,figsize=(9,3.4))
    markers={'dominance':'o','clearance':'s','lattice':'^','delay':'D','multidelay':'P'}
    for fam,marker in markers.items():
        oo=[r for r in oracle if r['family']==fam]
        x=[r['best']['panels'] for r in oo]
        y=[by[(r['name'],'local')]['panels'] for r in oo]
        axes[0].scatter(x,y,s=30,marker=marker,label=fam,alpha=.8)
    lim=max(max(r['best']['panels'] for r in oracle),max(r['panels'] for r in rows if r['family']!='shift'))
    axes[0].plot([4,lim],[4,lim],'k--',lw=.8)
    axes[0].set(xscale='log',yscale='log',xlabel='Best fixed center among nine: panels',ylabel='Local gauge: panels')
    axes[0].legend(fontsize=8,frameon=False)
    # Replay speed vs construction overhead relative to midrange.
    rt={(r['name'],r['policy']):r for r in replay}
    for fam,marker in markers.items():
        rr=[r for r in rows if r['family']==fam and r['policy']=='local']
        axes[1].scatter([by[(r['name'],'global')]['seconds']/r['seconds'] for r in rr],
            [rt[(r['name'],'global')]['seconds']/rt[(r['name'],'local')]['seconds'] for r in rr],
            s=25,marker=marker,alpha=.8)
    axes[1].axvline(1,color='k',ls='--',lw=.8);axes[1].axhline(1,color='k',ls='--',lw=.8)
    axes[1].set(xscale='log',yscale='log',xlabel='Construction speed ratio: midrange / local',ylabel='Replay speed ratio: midrange / local')
    fig.tight_layout();save(fig,'paired')
    # 4: finite stability atlas; exact counts checked independently by quadrature.
    fig,axes=plt.subplots(1,2,figsize=(9,3.4))
    taus=[.5,1,2,4];bs=[.5,2,4]
    arr=np.array([[next(r['count'] for r in rows if r['family']=='delay' and r['policy']=='local'
                  and float(Q(r['tau']))==t and float(Q(r['b']))==b) for t in taus] for b in bs])
    im=axes[0].imshow(arr,cmap='YlGnBu',vmin=0,vmax=6)
    for i in range(3):
        for j in range(4):axes[0].text(j,i,str(arr[i,j]),ha='center',va='center',color='white' if arr[i,j]>3 else 'black')
    axes[0].set(xticks=range(4),xticklabels=taus,yticks=range(3),yticklabels=bs,xlabel='Delay $\\tau$',ylabel='Feedback gain $b$',title='Certified unstable roots')
    cc=[r for r in rows if r['family']=='clearance' and r['policy']=='local']
    for p in ['raw','global','local']:
        rr=[r for r in rows if r['family']=='clearance' and r['policy']==p]
        axes[1].plot([int(r['name'].split('_')[1]) for r in rr],[r['panels'] for r in rr],'o-',color=COLORS[p],label=LABEL[p])
    axes[1].set(xscale='log',xlabel='Inverse boundary clearance $1/\\epsilon$',ylabel='Panels',title='Scaling does not remove nearby zeros')
    axes[1].legend(frameon=False);fig.tight_layout();save(fig,'application')
    # Tables and exact aggregate statistics.
    summary=[]
    for fam in markers:
        rr=[r for r in rows if r['family']==fam and r['policy']=='local']
        names=[r['name'] for r in rr]
        om={r['name']:r for r in oracle}
        result={'family':fam,'n':len(rr)}
        for p in ['raw','anchor','global','local']:
            result[p+'_median_panels']=statistics.median(by[(n,p)]['panels'] for n in names)
        result['oracle_median_panels']=statistics.median(om[n]['best']['panels'] for n in names)
        result['local_oracle_ratio']=statistics.median(om[n]['best']['panels']/by[(n,'local')]['panels'] for n in names)
        result['construction_ratio']=statistics.median(by[(n,'global')]['seconds']/by[(n,'local')]['seconds'] for n in names)
        result['replay_ratio']=statistics.median(rt[(n,'global')]['seconds']/rt[(n,'local')]['seconds'] for n in names)
        summary.append(result)
    nonshift=[r for r in rows if r['family']!='shift' and r['policy']=='local']
    om={r['name']:r for r in oracle}
    aggregate={'cases':len(rows)//4,'certificates':len(rows),
       'certified':sum(r['status']=='certified' for r in rows),
       'oracle_local_wins':sum(r['panels']<om[r['name']]['best']['panels'] for r in nonshift),
       'oracle_local_ties':sum(r['panels']==om[r['name']]['best']['panels'] for r in nonshift),
       'oracle_local_losses':sum(r['panels']>om[r['name']]['best']['panels'] for r in nonshift),
       'families':summary}
    (ROOT/'results/summary.json').write_text(json.dumps(aggregate,indent=2))
    lines=['\\begin{tabular}{lrrrrrr}\\toprule','Family & $n$ & None & Anchor & Midrange & Fixed-9 & Local\\\\\\midrule']
    for r in summary:
        lines.append(f"{r['family'].capitalize()} & {r['n']} & {r['raw_median_panels']:g} & {r['anchor_median_panels']:g} & {r['global_median_panels']:g} & {r['oracle_median_panels']:g} & {r['local_median_panels']:g}\\\\")
    lines+=['\\bottomrule\\end{tabular}']
    (ROOT/'paper/panels_table.tex').write_text('\n'.join(lines))
    lines=['\\begin{tabular}{lrrr}\\toprule','Family & Fixed-9 / local panels & Construction ratio & Replay ratio\\\\\\midrule']
    for r in summary:
        lines.append(f"{r['family'].capitalize()} & {r['local_oracle_ratio']:.2f} & {r['construction_ratio']:.2f} & {r['replay_ratio']:.2f}\\\\")
    lines+=['\\bottomrule\\end{tabular}']
    (ROOT/'paper/timing_table.tex').write_text('\n'.join(lines))
    print(json.dumps(aggregate,indent=2))

if __name__=='__main__':main()
