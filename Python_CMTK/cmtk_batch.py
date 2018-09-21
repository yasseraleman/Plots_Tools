import datetime
import time
import shutil
import os
import cmp, cmp.gui, cmp.connectome

cmpgui = cmp.gui.CMPGUI()

# DEFINE THE PROJECT
# ------------------
cmpgui.project_name = 'new_project_name'
cmpgui.project_dir = '/home/user/subjects/'
cmpgui.project_metadata = {}
print 'Project ' + cmpgui.project_name

# PARAMETERS
# ----------
# Actions
cmpgui.active_createfolder = True
cmpgui.active_dicomconverter = True
cmpgui.active_segmentation = True
cmpgui.active_registration = True
cmpgui.active_parcellation = True
cmpgui.active_applyregistration = True
cmpgui.active_reconstruction = True
cmpgui.active_tractography = True
cmpgui.active_fiberfilter = True
cmpgui.active_connectome = True
cmpgui.active_rsfmri = True
cmpgui.active_cffconverter = False
cmpgui.skip_completed_stages = False

# Configuration
cmpgui.freesurfer_home = '/usr/local/bin/freesurfer'
cmpgui.fsl_home = '/usr/share/fsl/4.1'
cmpgui.dtk_home = '/usr/local/bin/dtk'
cmpgui.fsloutputtype = 'NIFTI' # for fslview startup
cmpgui.pipeline_status_file = 'cmp.status'     # this file stores descriptions of the inputs/outputs to each stage of the CMP pipeline 

# Metadata
cmpgui.creator = u'Oscar Esteban'
cmpgui.description = ''
cmpgui.email = u'oesteban@die.upm.es'
cmpgui.emailnotify = []
cmpgui.license = ''
cmpgui.modified = None
cmpgui.publisher = u'BIT UPM'
cmpgui.reference = ''
cmpgui.relation = ''
cmpgui.rights = ''
cmpgui.species = 'Homo sapiens'

# DICOM converter
cmpgui.do_convert_T1 = True
cmpgui.do_convert_T2 = False
cmpgui.do_convert_diffusion = True
cmpgui.do_convert_fMRI = True
cmpgui.extract_diffusion_metadata = False
cmpgui.subject_raw_glob_T1 = '*.*'
cmpgui.subject_raw_glob_T2 = '*.*'
cmpgui.subject_raw_glob_diffusion = '*.*'

# Registration
cmpgui.registration_mode = 'BBregister' # 'Linear' / 'Nonlinear / BBregister'
#cmpgui.lin_reg_param = '-usesqform -nosearch -dof 6 -cost mutualinfo'
#cmpgui.nlin_reg_bet_T2_param = '-f 0.35 -g 0.15'
#cmpgui.nlin_reg_bet_b0_param = '-f 0.2 -g 0.2'
#cmpgui.nlin_reg_fnirt_param = '--subsamp=8,4,2,2 --mit... --applyrefmask=0,0,1,1'
cmpgui.bb_reg_param = '--init-header --dti'

# Segmentation 
cmpgui.recon_all_param = '-all -no-isrunning'
cmpgui.wm_handling = 1

# Parcellation
cmpgui.parcellation_scheme = 'NativeFreesurfer'

# Reconstruction
cmpgui.nr_of_b0 = str(1)

# Tractography
cmpgui.streamline_param = '--angle 60  --seeds 32'
cmpgui.streamline_param_dti = '--angle 60  --seeds 32'

# Fiber filtering
cmpgui.apply_fiberlength = True
cmpgui.apply_splinefilter = True
cmpgui.fiber_cutoff_lower = 8.0
cmpgui.fiber_cutoff_upper = 500.0

# Connectome creation
cmpgui.compute_curvature = False

# Measures in connectome
cmpgui.connection_P0 = True
cmpgui.connection_gfa = True
cmpgui.connection_kurtosis = True
cmpgui.connection_skewness = True
cmpgui.connection_adc = False
cmpgui.connection_fa = False

# Resting state - fMRI
cmpgui.rsfmri_registration_mode = 'BBregister' # 'Linear' / 'BBregister'
#cmpgui.rsfmri_lin_reg_param = '-usesqform -nosearch -dof 6 -cost mutualinfo'
cmpgui.rsfmri_bb_reg_param = '--init-header --dti'
cmpgui.do_save_mat = True;

# CFF converter
cmpgui.cff_cmatpickle = True
cmpgui.cff_fiberarr = True
cmpgui.cff_filteredfibers = True        
cmpgui.cff_finalfiberlabels = True       
cmpgui.cff_fullnetworkpickle = True      
cmpgui.cff_originalfibers = True        
cmpgui.cff_rawT1 = True                  
cmpgui.cff_rawT2 = True                  
cmpgui.cff_rawdiffusion = True           
cmpgui.cff_roisegmentation = True        
cmpgui.cff_scalars = True                
cmpgui.cff_surfacelabels = True          
cmpgui.cff_surfaces = True    

# LOOP THROUGH ALL THE SUBJECTS IN PROJECT DIRECTORY 
# --------------------------------------------------
for subj in os.listdir(cmpgui.project_dir):
	if os.path.isdir(os.path.join(cmpgui.project_dir,subj)):
		if subj.find('PH') != -1 and subj.find('PH0082') == -1 and subj.find('PH0080') == -1:	# SUBJECTS SELECTION!!!

			# Parameters for the current subject
			cmpgui.subject_name = subj
			cmpgui.subject_metadata = []
			print 'SUBJECT ' + cmpgui.subject_name

			# LOOP THROUGH ALL THE TIMEPOINTS OF THE CURRENT SUBJECTS
			# -------------------------------------------------------
			main_dir = os.path.join(cmpgui.project_dir, cmpgui.subject_name)

			for tp in os.listdir(main_dir):
				if os.path.isdir(os.path.join(main_dir,tp)) and (tp.find('scan1') != -1):	# TIME POINTS SELECTION

					# Parameters for the current time point
					cmpgui.subject_timepoint = tp
					cmpgui.subject_workingdir = os.path.join(main_dir, cmpgui.subject_timepoint)
					cmpgui.subject_logger = None # ?
					print 'TIME POINT ' + tp

					# Reconstruciton settings
					# DSI case
					cmpgui.diffusion_imaging_model = 'DSI'
					cmpgui.nr_of_sampling_directions = str(181)
					cmpgui.nr_of_gradient_directions = str(515)

#					# DTI case
#					elif tp.find('DTI') != -1:
#						cmpgui.diffusion_imaging_model = 'DTI'
#						cmpgui.max_b0_val = str(1000)
#						if tp.find('20') != -1:
#							cmpgui.gradient_table = 'siemens_20'
#						elif tp.find('64') != -1:
#							cmpgui.gradient_table = 'siemens_64'
#						else:
#							print cmpgui.subject_workingdir + ' does not follow siemens_20/64 scheme: SKIPPED!'
#							continue
#
#					# QBALL case
#					elif tp.find('QBI') != -1:
#						cmpgui.diffusion_imaging_model = 'QBALL'
#						cmpgui.gradient_table = 'siemens_256'
#						cmpgui.nr_of_gradient_directions = str(257)
#						cmpgui.nr_of_sampling_directions = str(181)

					print 'Processing ' + cmpgui.subject_name + ', ' + tp
					cmp.connectome.mapit(cmpgui)



