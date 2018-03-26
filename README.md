# Spinstractor
Scripts to extract and analyse spin evolution of black holes in Arepo snapshots

Files
-----

   - BH_Physics.py: master module with different functions and classes to work with simulated BHs as well as semi-analytical recipies for BH evolution.
   - BH_Indexing.py: script to detect unique BH IDs in simulations, regardless their redshift. **Produces:** *massdistro.txt, BH_IDs.txt*.
   - BH_Spin_Histories.py: script to extract BH spin, mass accretion and angular momentum histories from simulation for the BHs listed as unique. **Produces:** *spins/BH_(id).txt*.
   - BH_Mergers.py: script to extract descendant list for BH mergers from simulations. **Produces:** *BH_Mergers.txt*.
   
Format of output files
----------------------

   - massdistro.txt: all BHs at all times.
       - BH ID: BH id from simulation.
       - a(t): scale factor.
       - a_s: spin parameter.
       - BH mass: in units of 1e10 M_sun.
       - Spin mode: spin evolution mode (0 - prolonged, 1 - self-gravity).
    
   - BH_IDs.txt: BHs last appearance in simulation.
       - ID unique: unique ID assigned to BHs, instead of the long internal ID of the simulation. This is mass-ordered.
       - BH_ID: the long boring internal ID of the simulation.
       - a(t): scale factor at which this BH seen for the last time.
       - BH mass: in units of 1e10 M_sun.
       - a_s: spin parameter.
       - Spin mode: spin evolution mode (0 - prolonged, 1 - self-gravity).
       
   - spins/BH_(id).txt: spin evolution for the id-th BH.
       - a(t): scale factor.
       - a_s: spin parameter.
       - BH mass: in units of 1e10 M_sun.
       - dot_BH_mass: mass accretion rate in units of Eddington accretion.
       - spin_x: x-component of spin unit vector.
       - spin_y: y-component of spin unit vector.
       - spin_z: z-component of spin unit vector.
       - AM_x: x-component of angular momentum of gas around the BH.
       - AM_y: y-component of angular momentum of gas around the BH.
       - AM_z: z-component of angular momentum of gas around the BH.
       - Spin mode: spin evolution mode (0 - prolonged, 1 - self-gravity).
       
   - BH_Mergers.txt: list of descendants for each BH.
       - a(t): scale factor at which the merger happens.
       - BH 1 ID: ID_sim of BH main progenitor (remnant). Not necessarily more massive.
       - BH 1 mass: in units of 1e10 M_sun.
       - BH 2 ID: ID_sim of BH secondary progenitor (it disappears). Not necessarily more massive.
       - BH 2 mass: in units of 1e10 M_sun.
       
