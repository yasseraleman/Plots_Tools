function script_2_plot_multiple_surfaces_and_overlays
% aparc
SubjId = 'ch2';
hemi = 'lh';
Sfile = 'inflated';
cl = 'jet';
plotvar = 1;
CFile = '.aparc.annot'
OutDir = '/home/yaleman/picste'
hf = Plot_FreeS_results(SubjId,hemi,Sfile,CFile,cl);
Figurename = [OutDir filesep hemi '.' Sfile CFile '.jpg'];
a = getframe(hf);
imwrite(a.cdata,Figurename,'JPEG')
close(hf);

SubjId = 'ch2';
hemi = 'lh';
Sfile = 'sphere';
cl = 'jet';
plotvar = 1;
CFile = '.aparc.annot'
hf = Plot_FreeS_results(SubjId,hemi,Sfile,CFile,cl);
Figurename = [OutDir filesep hemi '.' Sfile CFile '.jpg'];
a = getframe(hf);
imwrite(a.cdata,Figurename,'JPEG')
close(hf);

SubjId = 'ch2';
hemi = 'lh';
Sfile = 'pial';
cl = 'jet';
plotvar = 1;
CFile = '.aparc.annot'
hf = Plot_FreeS_results(SubjId,hemi,Sfile,CFile,cl);
Figurename = [OutDir filesep hemi '.' Sfile CFile '.jpg'];
a = getframe(hf);
imwrite(a.cdata,Figurename,'JPEG')
close(hf);


% Thickness
SubjId = 'ch2';
hemi = 'lh';
Sfile = 'inflated';
cl = 'jet';
plotvar = 1;
CFile = '.thickness'
OutDir = '/home/yaleman/picste'
hf = Plot_FreeS_results(SubjId,hemi,Sfile,CFile,cl);
Figurename = [OutDir filesep hemi '.' Sfile CFile '.jpg'];
a = getframe(hf);
imwrite(a.cdata,Figurename,'JPEG')
close(hf);

SubjId = 'ch2';
hemi = 'lh';
Sfile = 'sphere';
cl = 'jet';
plotvar = 1;
CFile = '.thickness'
hf = Plot_FreeS_results(SubjId,hemi,Sfile,CFile,cl);
Figurename = [OutDir filesep hemi '.' Sfile CFile '.jpg'];
a = getframe(hf);
imwrite(a.cdata,Figurename,'JPEG')
close(hf);

SubjId = 'ch2';
hemi = 'lh';
Sfile = 'pial';
cl = 'jet';
plotvar = 1;
CFile = '.thickness'
hf = Plot_FreeS_results(SubjId,hemi,Sfile,CFile,cl);
Figurename = [OutDir filesep hemi '.' Sfile CFile '.jpg'];
a = getframe(hf);
imwrite(a.cdata,Figurename,'JPEG')
close(hf);


% Curvature
SubjId = 'ch2';
hemi = 'lh';
Sfile = 'inflated';
cl = 'hot';
plotvar = 1;
CFile = '.curv'
OutDir = '/home/yaleman/picste'
hf = Plot_FreeS_results(SubjId,hemi,Sfile,CFile,cl);
Figurename = [OutDir filesep hemi '.' Sfile CFile '.jpg'];
a = getframe(hf);
imwrite(a.cdata,Figurename,'JPEG')
close(hf);

SubjId = 'ch2';
hemi = 'lh';
Sfile = 'sphere';
cl = 'hot';
plotvar = 1;
CFile = '.curv'
hf = Plot_FreeS_results(SubjId,hemi,Sfile,CFile,cl);
Figurename = [OutDir filesep hemi '.' Sfile CFile '.jpg'];
a = getframe(hf);
imwrite(a.cdata,Figurename,'JPEG')
close(hf);

SubjId = 'ch2';
hemi = 'lh';
Sfile = 'pial';
cl = 'hot';
plotvar = 1;
CFile = '.curv'
hf = Plot_FreeS_results(SubjId,hemi,Sfile,CFile,cl);
Figurename = [OutDir filesep hemi '.' Sfile CFile '.jpg'];
a = getframe(hf);
imwrite(a.cdata,Figurename,'JPEG')
close(hf);

% local lgi
SubjId = 'ch2';
hemi = 'lh';
Sfile = 'inflated';
cl = 'jet';
plotvar = 1;
CFile = '.pial_lgi'
OutDir = '/home/yaleman/picste'
hf = Plot_FreeS_results(SubjId,hemi,Sfile,CFile,cl);
Figurename = [OutDir filesep hemi '.' Sfile CFile '.jpg'];
a = getframe(hf);
imwrite(a.cdata,Figurename,'JPEG')
close(hf);


SubjId = 'ch2';
hemi = 'lh';
Sfile = 'sphere';
cl = 'jet';
plotvar = 1;
CFile = '.pial_lgi'
OutDir = '/home/yaleman/picste'
hf = Plot_FreeS_results(SubjId,hemi,Sfile,CFile,cl);
Figurename = [OutDir filesep hemi '.' Sfile CFile '.jpg'];
a = getframe(hf);
imwrite(a.cdata,Figurename,'JPEG')
close(hf);



SubjId = 'ch2';
hemi = 'lh';
Sfile = 'pial';
cl = 'jet';
plotvar = 1;
CFile = '.pial_lgi'
hf = Plot_FreeS_results(SubjId,hemi,Sfile,CFile,cl);
Figurename = [OutDir filesep hemi '.' Sfile CFile '.jpg'];
a = getframe(hf);
imwrite(a.cdata,Figurename,'JPEG')
close(hf);
