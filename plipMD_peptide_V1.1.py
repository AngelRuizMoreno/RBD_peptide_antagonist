import MDAnalysis as md

import pandas as pd

from prody import confProDy,parsePDB,writePDB
silence_prody=confProDy(verbosity='none')

from plip.structure.preparation import PDBComplex

import progressbar

import os
import sys
import warnings
warnings.filterwarnings("ignore")


def plipmd(topol=None,traj=None):

	traj=list(traj.strip('[]').split(','))
	
	u = md.Universe(topol,traj)

	print ('\nINFO: your system contains {} segments with labels {} \n'.format(len(u.segments),list(u.segments.segids)))

	if len (u.segments.segids) ==1:
		print ('''
			WARNING: Only one segment was identified. Proceeding to manual segment definition\n\n

			You must define the peptide residues using MDAnalysis format. 
			
			For instance: resid 1:20
			
			This means peptide is from residue 1 to 20. Only peptide residues are required.

			Peptide residues will be added to "P" chain.

			Extra residues will be consider "RECEPTOR" (R). \n\n
			''')
		
		pep_residues=input('Type the peptide residues (example - resid 1:20 -):')
		pep_segment = u.add_Segment(segid='P')
		pep = u.select_atoms(pep_residues)
		pep.residues.segments=pep_segment

		rec_segment = u.add_Segment(segid='R')
		rec=u.select_atoms('protein and (not segid P)')
		rec.residues.segments=rec_segment

		print ('\nINFO: your system contains {} segments with labels {} \n'.format(len(u.segments),list(u.segments.segids)))

	
	peptide_segment=input ('\n\n1) Type the segment of your peptide:\n>')

	receptor_segment=input ('\n\n1) Type the segment of your receptor:\n>')

	sol_name=input ('\n2) Type the ResName of your water model. Must be 3-4 letter code (example - WAT or SOL or TIP3 -):\n>')

	for res in u.residues:
		if res.resname==sol_name:
			res.resname='HOH'
		if 'HI' in res.resname or 'HSD' in res.resname:
			res.resname='HIS'
		if 'CY' in res.resname:
			res.resname='CYS'
	for atom in u.atoms:
		if atom.name=='OH2':
			atom.name='OW'

	System=u.select_atoms('protein or resname HOH',updating=True)
	System=System.select_atoms('segid {} or (around 8 segid {})'.format(peptide_segment,peptide_segment),updating=True)

	table=pd.DataFrame()
	index=0
	print ('\nINFO: Your trajectory lenght is:{} steps\n'.format(range(len(u.trajectory))))
	start=int(input('4) Type the starting STEP to analyze:\n>'))
	finish=int(input('\n5) Type the ending STEP to analyze:\n>'))
	bar=progressbar.ProgressBar(max_value=finish)
	print ('\n\n-----  -----  -----  RUNNING THE ANALYSIS  -----  -----  -----\n\n')
	
	for i in range(start,finish):
		try:
			name='frame_tmp.pdb'
			PDB= md.Writer(name, multiframe=False)
			for ts in u.trajectory[i:i+1]:
				PDB.write(System)

				#Prody for converting peptide ATOM to HETATM
				pdb_file = parsePDB(name)
				labels=[True if i.getSegname()==peptide_segment else False for i in pdb_file]
				pdb_file.setFlags('hetatm',labels)
				writePDB(name,pdb_file)
				
				plip_job = PDBComplex()
				plip_job.load_pdb(name) 
				plip_job.analyze()
				
				for key in plip_job.interaction_sets:
					if key.split(':')[1]==peptide_segment:
						for interaction in plip_job.interaction_sets.get(key).all_itypes:
							if interaction.reschain != peptide_segment:
								table.loc[index,'Frame']=ts.frame
								table.loc[index,'Time']=ts.time
								interaction_type=str(type(interaction)).split('.')[-1].replace("'>","")
								table.loc[index,'Type']=interaction_type
								table.loc[index,'Protein']=interaction.restype+str(interaction.resnr)+'_'+interaction.reschain
								table.loc[index,'Peptide']=interaction.restype_l+str(interaction.resnr_l)+'_'+interaction.reschain_l
								
								if interaction_type == 'hbond':
									table.loc[index,'Acceptor']=interaction.atype
									table.loc[index,'AcceptorIdx']=interaction.a.idx
									table.loc[index,'Donor']=interaction.dtype
									table.loc[index,'DonorIdx']=interaction.d.idx
									table.loc[index,'DistanceAD']=interaction.distance_ad
									table.loc[index,'DistanceAH']=interaction.distance_ah
									table.loc[index,'Angle']=interaction.angle
									table.loc[index,'Force']=interaction.type
									table.loc[index,'ProtIsDon']=interaction.protisdon
								
								elif interaction_type == 'pication':
									table.loc[index,'Charge']=interaction.charge.type
									table.loc[index,'ChargedAtoms']=",".join([i.type for i in interaction.charge.atoms])
									table.loc[index,'Force']=interaction.type
									table.loc[index,'RingType']=interaction.ring.type
									table.loc[index,'RingAtoms']=",".join([i.type for i in interaction.ring.atoms])
									table.loc[index,'RingAtomsIdx']=",".join([str(i.idx) for i in interaction.ring.atoms])        
								
								elif interaction_type=='saltbridge':
									table.loc[index,'NegAtoms']=",".join([i.type for i in interaction.negative.atoms])
									table.loc[index,'NegAtomsIdx']=",".join([str(i.idx) for i in interaction.negative.atoms])
									table.loc[index,'PosAtoms']=",".join([i.type for i in interaction.positive.atoms])
									table.loc[index,'PosAtomsIdx']=",".join([str(i.idx) for i in interaction.positive.atoms])
									table.loc[index,'Distance']=interaction.distance
									table.loc[index,'ProtIsPos']=interaction.protispos
									
								elif interaction_type == 'hydroph_interaction':
									table.loc[index,'RecAtom']=interaction.bsatom.type
									table.loc[index,'RecAtomIdx']=interaction.bsatom.idx
									table.loc[index,'LigAtom']=interaction.ligatom.type
									table.loc[index,'LigAtomIdx']=interaction.ligatom.idx
									table.loc[index,'Distance']=interaction.distance
									
								elif interaction_type == 'pistack':
									table.loc[index,'RecRingType']=interaction.proteinring.type
									table.loc[index,'LigRingType']=interaction.ligandring.type
									table.loc[index,'RecRingAtoms']=",".join([i.type for i in interaction.proteinring.atoms])
									table.loc[index,'RecAtomsIdx']=",".join([str(i.idx) for i in interaction.proteinring.atoms])
									table.loc[index,'LigRingAtoms']=",".join([i.type for i in interaction.ligandring.atoms])
									table.loc[index,'LigRingAtomsIdx']=",".join([str(i.idx) for i in interaction.ligandring.atoms])
									table.loc[index,'Distance']=interaction.distance
									table.loc[index,'Angle']=interaction.angle
									table.loc[index,'Offset']=interaction.offset
									
								elif interaction_type == 'waterbridge':
									table.loc[index,'AccType']=interaction.atype
									table.loc[index,'DonType']=interaction.dtype
									table.loc[index,'WaterIdx']=interaction.water_orig_idx
									table.loc[index,'DistanceAWat']=interaction.distance_aw
									table.loc[index,'DistanceDWat']=interaction.distance_dw
									table.loc[index,'AngleDon']=interaction.d_angle
									table.loc[index,'AngleWat']=interaction.w_angle
									table.loc[index,'ProtIsDon']=interaction.protisdon

								elif interaction_type == 'halogenbond':
									table.loc[index,'Acceptor']=interaction.acctype
									table.loc[index,'Donor']=interaction.acctype
									table.loc[index,'Distance']=interaction.distance
									table.loc[index,'DonAngle']=interaction.don_angle
									table.loc[index,'AccAngle']=interaction.acc_angle
					  
								elif interaction_type=='metal_complex':
									table.loc[index,'MetalType']=interaction.metal.type
									table.loc[index,'Idx']=interaction.metal.idx
									table.loc[index,'TargetType']=interaction.target_type
									table.loc[index,'FunctGroup']=interaction.target.fgroup
									table.loc[index,'Geometry']=interaction.geometry
									table.loc[index,'Distance']=interaction.distance
									table.loc[index,'Location']=interaction.location
								
								index=index+1    
			bar.update(i+1)
			os.remove(name)
		
		except Exception:
			continue

	print ('\n\n-----  -----  -----  SAVING THE RESULTS, PLEASE WAIT  -----  -----  -----\n\n')	
	table.set_index(['Frame','Time'], inplace=True)
	table.sort_index(inplace=True)
	table.to_excel('Interactions_Table.xlsx')
	print ('\n\n***** ***** ***** ALL DONE, DATA SAVED ON: Interactions_Table.xlsx ***** ***** *****\n\n')	 
if __name__ == "__main__":
	
	plipmd(sys.argv[1],sys.argv[2])
