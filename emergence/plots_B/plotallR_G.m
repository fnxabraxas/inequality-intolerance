
function plotallR_G(i,stgv,ib,bnam,icont)

  %dirs=strvcat( 'broteR_TT', 'broteR_TR', 'broteT_RT', 'broteT_RR'); %,'broteT_RT_PD');
  %%%dirs=strvcat( 'broteR_TT'); %, 'broteR_TR', 'broteT_RT', 'broteT_RR');
  dirlon=[9 9 9 9];
  %colordir=strvcat( 'b', 'y', 'r', 'r' );
  %%%colordir=[ 0 0.3 1; 1 1 0; 0.8 0.03 0.03; 0.8 0.03 0.03 ];
  %%%transpdir=[1 0.5 1.0 0.5];
  lindir= strvcat( 'none', 'none', 'none', 'r' );
  %%%cambiaY=[0 1 0 0];
  figure(i);

  dirs=strvcat( 'broteR_TT', 'broteT_RT', 'broteR_TR' ); %, 'broteT_RR' )
  colordir=[ 0 0.3 1; 0.8 0.03 0.03; 1 1 0; 0 0 0 ];
  transpdir=[1 0.5 0.5 0.5];
  cambiaY=[0 0 1 0];
  
  %if ib>=3
  %  ib2=ib+4
  %else
  %  ib2=ib*2-1
  %end
  %subplot(2,5,[ib2 ib2+1]);
  
  switch icont
    case 1,
      subplot('Position',[0.1 0.66 0.22 0.29]);
    case 2,
      subplot('Position',[0.415 0.66 0.22 0.29]);
    case 3,
      subplot('Position',[0.73 0.66 0.22 0.29]);
    case 4,
      subplot('Position',[0.1 0.23 0.22 0.29]);
    case 5,
      subplot('Position',[0.415 0.23 0.22 0.29]);
  end
  

  %axis(axis)
  axis ([-1 1 0 1]);
  
  ylabel('1-\epsilon');
  %ylabel('G_A');
  xlabel('\epsilon_B-\epsilon_A');
  hold all;

  [Ndir,tt]=size(dirs);

  %Ndir=1
  for idir=1:1 %Ndir
    %file=['../modif/' dirs(idir,1:dirlon(idir)) '_101_' num2str(stgv(i)) '.dat'];
    if(cambiaY(idir)==1)
      file=['../modif_B/' dirs(idir,1:dirlon(idir)) '_' bnam '_' num2str(stgv(i)) '_rev_101.dat']  
    else
      file=['../modif_B/' dirs(idir,1:dirlon(idir)) '_' bnam '_' num2str(stgv(i)) '_101.dat']
    end
    disp(file);
    Cfile=importdata(file,' ',1);
    [Ny,Nb]=size(Cfile.data);
    Ny=Ny-1;
    Nb=(Nb-1)/2;
    epsB=Cfile.data(2:Ny+1,1);
    epsA=Cfile.data(2:Ny+1,ib+1);
    epsBI=Cfile.data(2:Ny+1,ib+Nb+1);
    y=Cfile.data(1,ib+1);
    title (['y=' num2str(y)]);
    %if stgv(i)~=405 && stgv(i)~=469  plot([-1 2],[1-1/bhead 1-1/bhead],'b-.'); end
    
    xm=[-1:0.01:1];
    ym1A=1-xm*(1-y);
    ym2A=-xm*(1-y);
    ym1B=1+xm*y;
    ym2B=xm*y;
    plot(xm,ym1A,'k-');
    plot(xm,ym2A,'k-');
    plot(xm,ym1B,'k-');
    plot(xm,ym2B,'k-');
    
    fill([-1 -1 0],[0 1-y 0],[0.5 0.5 0.5]);
    fill([-1 -1 0],[1-y 1 1],[0.5 0.5 0.5]);
    fill([1 1 0],[y 0 0],[0.5 0.5 0.5]);
    fill([1 1 0],[1 y 1],[0.5 0.5 0.5]);
    
    epsTt=y*epsA+(1-y)*epsB;
    GAt=(epsB-epsA);
    %GAt=(1-y)*epsB-y*epsA;
    [epsT,id]=sort(epsTt); GA=GAt(id);
    plot(GA,1-epsT,'b-');
    %plot([0:0.01:1],GAmax,'k--'); plot([0:0.01:1],GAmin,'k--');

  end

  %plot([0 0], [0 1], 'k:');
  plot([0 0], [0 1], 'k:');
      
  %axis ([-1 1 0.01 0.99]);

  box

end

