
%stgv=[ 405 469 1437 1501 2186 2202 2219 2234 2442 2453 2458 2475 2490 2517 3477 3541];


%stgv=[ 405 469 ];
%stgv=[ 1437 1501 2453 2517 3477 3541];




yallow=[2 4 6];


Nbb=length(yallow);
stgv=[ 1437 2517] %[ 1437 1501 2453 2517 3477 3541];
bnamevec={'2'} %{'1.5', '2', '3'};

for ibb=1:length(bnamevec)
for i=1:length(stgv)
 h=figure(i);
 set(h,'Position',[0,0,750,550]);
 clf reset;
    %break
 icont=0;
 for ib=1:Nbb
    icont=icont+1;
    
    if ((stgv(i)==1437)&&(bnamevec{ibb}=='2'))
      interv=[[-99 -99]; [0.817 0.48]; [-99 -99]; [0.514 0.33]; [0 0]; [0 0]; [0 0]; [0 0]; [0 0]; [0 0]; [0 0]];
    elseif ((stgv(i)==2517)&&(bnamevec{ibb}=='2'))
      interv=[[1 0]; [1 0]; [1 0]; [1 0]; [0 0]; [0 0]; [0 0]; [0 0]; [0 0]; [0 0]; [0 0]];    
    end
    
    plotallR_bi(i,stgv,yallow(ib),bnamevec{ibb},icont,interv(yallow(ib),:));
    %break
 end

  %plotLEG([1 1 1 0 1 0]);

  filename=['racism_bcte_' bnamevec{ibb} '_' num2str(stgv(i)) '_bi'];

  filenamesvg=[filename '.svg'];
  plot2svg(filenamesvg)


  filenameps=[filename '.eps'];
  
  %filenamepdf=[filename '.pdf'];
  print ('-dpsc2',filenameps);
  %fileps2pdf=['ps2pdf -sPAPERSIZE=a4 ' filenameps ' ' filenamepdf];
  %system(fileps2pdf);
  %zippdf=['zip racism_pdf ' filenamepdf];
  %system(zippdf);

  %%%filenameps=[filename '.pdf'];
  %%%print ('-dpdf','-r600',filenameps);

  %%%filenameps=[filename '.jpg'];
  %%%print ('-djpeg',filenameps);

  %%%filenameps=[filename '.bmp'];
  %%%print ('-dbmp',filenameps);

end
end

