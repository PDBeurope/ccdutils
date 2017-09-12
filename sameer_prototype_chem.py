# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

#from __future__ import print_function
import os
import re

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw.MolDrawing import MolDrawing
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Geometry import rdGeometry


from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
#from IPython.display import SVG

from mmCif import *
import mmCif.mmcifIO as mmcif
from utilities import supply_list_of_sample_cifs, fragment_library_file_path


#from Chem import rdchem
BondType = {  'SING': Chem.rdchem.BondType.SINGLE, 'DOUB': Chem.rdchem.BondType.DOUBLE,'AROM': Chem.rdchem.BondType.AROMATIC }
StereoType = { 'S': Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW, 'R':Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW}
ElementType = {'HE':'He','LI':'Li','BE':'Be','NE':'Ne','NA':'Na','FE': 'Fe'}

this_script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(this_script_dir, 'data')

if __name__ == "__main__":

#    mol1 = Chem.MolFromSmiles('c1ccccc1O')
    smiles={}
    FragmentMols = {}
    
    fragmentFile = fragment_library_file_path
    open(fragmentFile, 'r')
    
    lines = [line.strip("\n").strip("\t").strip("\r").replace("\n","").replace("\t","").replace("\r","") for line in open(fragmentFile)]
    smiles=dict(((smile, name) for name, smile in (line.split(':') for line in lines)))
#    print smiles
    for k in smiles:
        m = Chem.MolFromSmiles(k)
#        sanFragMol = Chem.SanitizeMol(m,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        FragmentMols[k]= m
    
    
#    print FragmentMols
    
    
    cif_parser = mmcif.CifFileReader(input='data', preserve_order=True)
   
    cif_in_dir = os.path.join(this_script_dir, 'tests', 'ccd_mmcif_test_files')
    sdf_out_dir = os.path.join(this_script_dir, 'tests', 'out', 'sdf_out')
    svg_out_dir = os.path.join(this_script_dir, 'tests', 'out', 'svg_out')
    for out_dir in [sdf_out_dir, svg_out_dir]:  # make output directories if necessary
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)
    imgFileName = ''
    
    for fileName in supply_list_of_sample_cifs():
        skip = ('00O', '03R', '0OD', '10R', '3CD', 'ASX', 'CDL', 'CMO', 'D3O', 'NA', 'UNL')
        if any(substring in fileName for substring in skip):
            continue
        # print('fileName={}'.format(fileName))

        mol2 = Chem.MolFromSmiles('')
        mol2.UpdatePropertyCache(strict=False)
        
        em = Chem.RWMol(mol2)
        PDBconf = Chem.Conformer()
        PDBconf.SetId(0)
        Idealconf = Chem.Conformer()
        Idealconf.SetId(1)
        cifObj = cif_parser.read(os.path.join(cif_in_dir, fileName), output='cif_wrapper')
        chem_comp = list(cifObj.values())[0]
        chem_comp_id = chem_comp._chem_comp['id'][0]
        print(chem_comp_id)
#        imgFileName = path2+chem_comp._chem_comp['id'][0]+".png"
#        print imgFileName
        atoms=chem_comp._chem_comp_atom
        bonds=chem_comp._chem_comp_bond 
        rdatoms=[]
        rdatomdict={}
#        aromatics = []
#        inRings = []
        for atom in atoms:
            atomName = atom['atom_id']
            atomType = atom['type_symbol']
            if atomType in ElementType:
                atomType = ElementType[atomType]
            

            id = atom['pdbx_ordinal']
            arom = atom['pdbx_aromatic_flag']
            ring = atom['pdbx_aromatic_flag']
            stereo = atom['pdbx_stereo_config']
            pdb_x = float(atom['model_Cartn_x'])
            pdb_y = float(atom['model_Cartn_y'])
            pdb_z = float(atom['model_Cartn_z'])
            # print(pdb_x," ",pdb_y," ",pdb_z)

            ideal_x = float(atom['pdbx_model_Cartn_x_ideal'])
            ideal_y = float(atom['pdbx_model_Cartn_y_ideal'])
            ideal_z = float(atom['pdbx_model_Cartn_z_ideal'])
            # print(ideal_x," ",ideal_y," ",ideal_z)



            
            rdAtom = Chem.Atom(atomType)
            rdAtom.SetProp('name',atomName)
            rdAtom.SetProp('ordinal',id)
#            rdAtom.SetProp('_CIPCode',stereo)
            
            if (stereo != 'N'):
                # print(StereoType[stereo])
                Chem.Atom.SetChiralTag(rdAtom,StereoType[stereo])
            
            
            rdPoint_PDB = rdGeometry.Point3D(pdb_x,pdb_y,pdb_z)
            rdPoint_Ideal = rdGeometry.Point3D(ideal_x,ideal_y,ideal_z)
            
            idx = em.AddAtom(rdAtom)
#            print idx, atomName
            PDBconf.SetAtomPosition(idx,rdPoint_PDB)
            Idealconf.SetAtomPosition(idx,rdPoint_Ideal)
            
            rdatoms.append(rdAtom)
            rdatomdict[atomName]= rdAtom
            

        for k in rdatoms:
            for searchResult in bonds.searchiter('atom_id_1', k.GetProp('name')):
#                print k.GetProp('ordinal')
#                print searchResult['atom_id_2']
#                print rdatomdict[searchResult['atom_id_2']].GetProp('ordinal')
#                print searchResult['value_order']
                aromaticBond = searchResult['pdbx_aromatic_flag']
#                print aromaticBond
                if aromaticBond == 'Y':
                    bondType = Chem.rdchem.BondType(BondType['AROM'])
                else:
                    bondType = Chem.rdchem.BondType(BondType[searchResult['value_order']])
                rdBond = em.AddBond(k.GetIntProp('ordinal')-1,rdatomdict[searchResult['atom_id_2']].GetIntProp('ordinal')-1,bondType)
#                print searchResult
#            print "******************"
#        m = em.GetRingInfo()
#        print m
#        
#        for atom in aromatics:
#            atom.SetIsAromatic(True)
        mol2 = em.GetMol()
        
            
       
        
  
  
        mol2.AddConformer(PDBconf)
        Chem.rdmolops.WedgeMolBonds(mol2,PDBconf)
        w_PDB = Chem.SDWriter(sdf_out_dir + '/' + chem_comp_id + '_PDB.sdf')
        w_PDB.SetKekulize(False)
        w_PDB.write(mol2, confId=0)
        
        mol2.AddConformer(Idealconf)
        Chem.rdmolops.WedgeMolBonds(mol2,Idealconf)
        w_Ideal = Chem.SDWriter(sdf_out_dir + '/' + chem_comp._chem_comp['id'][0] + '_Ideal.sdf')
        w_Ideal.SetKekulize(False)
        w_Ideal.write(mol2, confId=1)
#        print Chem.MolToMolBlock(mol2)
#        print Chem.SanitizeMol(mol2,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
#        print mol2
#        print Chem.MolToSmiles(mol2)
#        Draw.MolToFile(mol2,imgFileName)
        drawer_nH = rdMolDraw2D.MolDraw2DSVG(800,800)
        
#        drawer_nH.setFontSize(0.3)
        opts1 = drawer_nH.drawOptions()
#        Chem.SanitizeMol(mol2)



        
        
        for i in range(mol2.GetNumAtoms()):
#            print mol2.GetAtomWithIdx(i).GetProp('name')
            opts1.atomLabels[i] = mol2.GetAtomWithIdx(i).GetSymbol()
            
        confIdFor2D= AllChem.Compute2DCoords(mol2,canonOrient=True)
        
#        confIdFor2D= AllChem.GenerateDepictionMatching3DStructure(mol2, mol2, confId=0)
        drawer_nH.DrawMolecule(mol2,confId=confIdFor2D)
        drawer_nH.FinishDrawing()
#        svg = drawer.GetDrawingText().replace('svg:','')
        svgFileName = svg_out_dir + '/' + chem_comp._chem_comp['id'][0] + '_NoLable_H.svg'
        f = open(svgFileName,'w')
        f.write(drawer_nH.GetDrawingText().replace('svg:','').replace('xmlns:svg','xmlns'))



        drawer = rdMolDraw2D.MolDraw2DSVG(800,800)
        
        opts = drawer.drawOptions()
#        Chem.SanitizeMol(mol2)



        
        
        for i in range(mol2.GetNumAtoms()):
            opts.atomLabels[i] = mol2.GetAtomWithIdx(i).GetProp('name')

        confIdFor2D= AllChem.Compute2DCoords(mol2,canonOrient=True)
#        confIdFor2D= AllChem.GenerateDepictionMatching3DStructure(mol2, mol2, confId=0)
        drawer.DrawMolecule(mol2,confId=confIdFor2D)
        drawer.FinishDrawing()
#        svg = drawer.GetDrawingText().replace('svg:','')
        svgFileName = svg_out_dir + '/' + chem_comp._chem_comp['id'][0] + '_atomLable_H.svg'
        f = open(svgFileName,'w')
        f.write(drawer.GetDrawingText().replace('svg:','').replace('xmlns:svg','xmlns'))





        matches={} 
#        print em
        for frag in FragmentMols:
#            print frag
            if mol2.HasSubstructMatch(FragmentMols[frag]):
#                print FragmentMols[frag]
                for allIndex in mol2.GetSubstructMatches(FragmentMols[frag]):
                    atomname =[]
                    for index in allIndex:
                        atomname.append(mol2.GetAtomWithIdx(index).GetProp('name'))
#                    print '___________________'
                matches[smiles[frag]]=atomname
        
        print(matches)
        print('***************************')
    
    
