import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
from pathlib import Path

OUT = Path(os.environ.get("DOR_FIGURE_OUTPUT_DIR", ".")).resolve()
OUT.mkdir(parents=True, exist_ok=True)

def save_all(fig, stem, dpi=600):
    fig.savefig(OUT/f'{stem}.png', dpi=dpi, bbox_inches='tight')
    fig.savefig(OUT/f'{stem}.pdf', bbox_inches='tight')

# ---------------- Figure 1 ----------------
fig=plt.figure(figsize=(10,10))
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')

def box(x,y,w,h,title,lines, lw=1.1):
    p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.012,rounding_size=0.012',linewidth=lw,fill=False)
    ax.add_patch(p)
    ax.text(x+0.02*w,y+h-0.22*h,title,fontsize=12,fontweight='bold',va='center')
    for i,line in enumerate(lines):
        ax.text(x+0.02*w,y+h-(0.48+i*0.15)*h,line,fontsize=10,va='center')
    return p

def arrow(x1,y1,x2,y2):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='->',mutation_scale=9,linewidth=0.9))

ax.text(0.04,0.965,'1 COHORT INPUTS AND PROVENANCE',fontsize=14,fontweight='bold',va='top')
box(0.04,0.805,0.27,0.12,'GSE274832',['147 clinical context','6 omics / independent','3 DOR + 3 NOR','metadata conflict disclosed'])
box(0.365,0.805,0.27,0.12,'GSE193136',['24 recruited','12 omics / independent','6 DOR + 6 NOR','age-adjusted primary model'])
box(0.69,0.805,0.27,0.12,'GSE232306',['60 recruited','12 omics / independent','6 DOR + 6 NOR','covariates unavailable'])
for x in [0.175,0.50,0.825]: arrow(x,0.805,x,0.765)
core=FancyBboxPatch((0.16,0.73),0.68,0.045,boxstyle='round,pad=0.01',linewidth=1.0,fill=False)
ax.add_patch(core)
ax.text(0.50,0.753,'Core RNA-seq analysis: 30 independent transcriptomes\n15 DOR + 15 NOR; cohorts modeled separately',fontsize=10.5,fontweight='bold',ha='center',va='center')
box(0.04,0.65,0.40,0.055,'E-MTAB-391 - legacy sensitivity only',['28 cycle-samples / 26 subjects; not pooled into core synthesis'])
arrow(0.50,0.73,0.50,0.61)
ax.text(0.04,0.625,'2 COHORT-SPECIFIC EXECUTION',fontsize=14,fontweight='bold',va='top')
box(0.08,0.51,0.84,0.075,'Repeated independently for each core accession',['Species/tissue identity and cohort QC','Salmon/tximport quantification with cohort-specific DESeq2 models','DOR-minus-NOR effect tables plus leave-one-sample-out influence audits'])
ax.text(0.04,0.455,'3 SEQUENTIAL CROSS-COHORT SYNTHESIS',fontsize=14,fontweight='bold',va='top')
arrow(0.50,0.51,0.50,0.415); arrow(0.50,0.415,0.17,0.415); arrow(0.17,0.415,0.17,0.37)
box(0.04,0.28,0.27,0.09,'Gene-level effects',['Pairwise concordance','Heterogeneity assessment'])
box(0.365,0.28,0.27,0.09,'Pathway convergence',['Full-rank Hallmark','and Reactome scores','Prespecified convergence rule'])
box(0.69,0.28,0.27,0.09,'Held-out LOCO',['Retained-pair selection','Held-out strict criterion'])
arrow(0.31,0.325,0.365,0.325); arrow(0.635,0.325,0.69,0.325)
ax.text(0.04,0.235,'4 INTEGRATED INTERPRETATION',fontsize=14,fontweight='bold',va='top')
arrow(0.825,0.28,0.825,0.195); arrow(0.825,0.195,0.50,0.195); arrow(0.50,0.195,0.50,0.17)
box(0.08,0.105,0.84,0.075,'Integrated interpretation',['Cohort-dependent effects with weak gene-level reproducibility','Selective pathway convergence and a narrow internally stable pathway signal'])
arrow(0.50,0.105,0.50,0.085)
box(0.08,0.020,0.84,0.060,'Interpretation boundary',['Internal conditional stability across the same public cohorts','No biomarker, causal mechanism, external validation, or clinical claim'])
save_all(fig,'Figure_1')
plt.close(fig)

# ---------------- Figure 5 ----------------
fig=plt.figure(figsize=(14,11))
gs=fig.add_gridspec(2,2,hspace=0.62,wspace=0.62)
colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
rot=['Hold out\n274832','Hold out\n193136','Hold out\n232306']
width=0.22
pos_h=np.array([0.0,1.0,2.0]); pos_r=np.array([4.1,5.1,6.1])

# A: counts
axA=fig.add_subplot(gs[0,0])
h_pair=[1,4,14]; h_dir=[1,4,5]; h_strict=[1,1,1]
r_pair=[30,25,36]; r_dir=[28,22,11]; r_strict=[7,7,7]
for off, vals_h, vals_r, lab, c in [
    (-width,h_pair,r_pair,'Pair-selected',colors[0]),
    (0,h_dir,r_dir,'Held-out direction',colors[1]),
    (width,h_strict,r_strict,'Held-out strict criterion',colors[2])]:
    axA.bar(pos_h+off,vals_h,width=width,label=lab,color=c)
    axA.bar(pos_r+off,vals_r,width=width,color=c)
    for p,v in zip(pos_h+off,vals_h): axA.text(p,v+0.7,str(v),ha='center',fontsize=8)
    for p,v in zip(pos_r+off,vals_r): axA.text(p,v+0.7,str(v),ha='center',fontsize=8)
axA.set_ylim(0,45)
axA.set_xticks(list(pos_h)+list(pos_r),rot+rot,fontsize=8)
axA.set_ylabel('Pathways')
axA.set_title('A  LOCO selection and held-out concordance',loc='left',fontweight='bold',pad=10)
axA.text(1.0,17.2,'HALLMARK',ha='center',fontweight='bold',fontsize=9)
axA.text(5.1,40.5,'REACTOME',ha='center',fontweight='bold',fontsize=9)
axA.legend(frameon=False,fontsize=8,loc='upper center',bbox_to_anchor=(0.5,-0.24),ncol=3)
axA.spines[['top','right']].set_visible(False)

# B: rates
axB=fig.add_subplot(gs[0,1])
h_direction=[100,100,36]; h_strict_pct=[100,25,7]
r_direction=[93,88,31]; r_strict_pct=[23,28,19]
for off, vals_h, vals_r, lab, c in [
    (-width/2,h_direction,r_direction,'Same direction',colors[0]),
    (width/2,h_strict_pct,r_strict_pct,'Strict criterion',colors[1])]:
    axB.bar(pos_h+off,vals_h,width=width,label=lab,color=c)
    axB.bar(pos_r+off,vals_r,width=width,color=c)
    for p,v in zip(pos_h+off,vals_h):
        dy = 8 if (lab == 'Strict criterion' and v >= 95) else 3
        axB.text(p,v+dy,f'{v}%',ha='center',fontsize=8)
    for p,v in zip(pos_r+off,vals_r): axB.text(p,v+3,f'{v}%',ha='center',fontsize=8)
axB.set_ylim(0,126)
axB.set_yticks([0,25,50,75,100],['0%','25%','50%','75%','100%'])
axB.set_xticks(list(pos_h)+list(pos_r),rot+rot,fontsize=8)
axB.set_ylabel('Rate among pair-selected pathways')
axB.set_title('B  Held-out criterion rates',loc='left',fontweight='bold',pad=10)
axB.text(1.0,116,'HALLMARK',ha='center',fontweight='bold',fontsize=9)
axB.text(5.1,116,'REACTOME',ha='center',fontweight='bold',fontsize=9)
axB.legend(frameon=False,fontsize=8,loc='upper center',bbox_to_anchor=(0.5,-0.24),ncol=2)
axB.spines[['top','right']].set_visible(False)

# C: NES heatmap
axC=fig.add_subplot(gs[1,0])
paths=['H: P53 Pathway','R: AUF1/HNRNPD destabilizes mRNA','R: Translation initiation','R: Nervous-system development','R: Neutrophil degranulation','R: RAS regulation by GAPs','R: rRNA processing','R: RUNX2 transcriptional regulation']
mat=np.array([[-1.80,-1.62,-1.58],[-2.02,-1.71,-2.29],[-1.49,-2.44,-3.34],[-1.46,-1.89,-1.79],[-1.75,-1.47,-1.50],[-1.69,-1.61,-1.65],[-2.22,-1.69,-3.00],[-1.63,-1.54,-1.54]])
im=axC.imshow(mat,aspect='auto')
axC.set_xticks([0,1,2],['274832','193136','232306'])
axC.set_yticks(range(len(paths)),paths,fontsize=8)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]): axC.text(j,i,f'{mat[i,j]:.2f}',ha='center',va='center',fontsize=8)
axC.set_title('C  Internally stable LOCO pathway set',loc='left',fontweight='bold',pad=10)
# compact in-panel NES scale, away from panel D
cbar=fig.colorbar(im,ax=axC,orientation='horizontal',fraction=0.055,pad=0.12,aspect=28)
cbar.set_label('NES',fontsize=8); cbar.ax.tick_params(labelsize=7)

# D: Hallmark pass matrix
axD=fig.add_subplot(gs[1,1])
h_paths=['P53 Pathway','Apoptosis','Coagulation','MYC Targets V1','MYC Targets V2','Myogenesis','ROS Pathway','UV Response Up']
passmat=np.array([[1,1,1],[0,0,1],[0,0,1],[0,1,0],[0,1,0],[0,0,1],[0,1,0],[0,0,1]])
axD.imshow(passmat,aspect='auto',vmin=0,vmax=1)
axD.set_xticks([0,1,2],['Without\n274832','Without\n193136','Without\n232306'],fontsize=8)
axD.set_yticks(range(len(h_paths)),h_paths,fontsize=8)
for i in range(passmat.shape[0]):
    for j in range(passmat.shape[1]): axD.text(j,i,'PASS' if passmat[i,j] else '-',ha='center',va='center',fontsize=8)
axD.set_title('D  Hallmark stability across LOCO rotations',loc='left',fontweight='bold',pad=10)

save_all(fig,'Figure_5')
plt.close(fig)
print('generated Figure_1 and Figure_5')
