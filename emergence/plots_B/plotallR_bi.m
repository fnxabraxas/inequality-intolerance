
function plotallR_bi(i,stgv,ib,bnam,icont,inter)


  %dirs=strvcat( 'broteR_TT', 'broteR_TR', 'broteT_RT', 'broteT_RR'); %,'broteT_RT_PD');
  %%%dirs=strvcat( 'broteR_TT'); %, 'broteR_TR', 'broteT_RT', 'broteT_RR');
  dirlon=[9 9 9 9];
  %colordir=strvcat( 'b', 'y', 'r', 'r' );
  %%%colordir=[ 0 0.3 1; 1 1 0; 0.8 0.03 0.03; 0.8 0.03 0.03 ];
  %%%transpdir=[1 0.5 1.0 0.5];
  lindir= strvcat( 'none', 'none', 'none', 'r' );
  %%%cambiaY=[0 1 0 0];
  figure(i);

  dirs=strvcat( 'broteR_TT') %, 'broteT_RT' )
  colordir=[ 0 0.3 1; 0.8 0.03 0.03; 1 1 0; 0.8 0.03 0.03 ];
  transpdir=[1 0.5 1.0 0.5];
  cambiaY=[0 0 ];

  xmin=0.01; xmax=0.99; ymin=0.01; ymax=0.99;
  
  %if ib>=3
  %  ib2=ib+4
  %else
  %  ib2=ib*2-1
  %end
  %subplot(2,5,[ib2 ib2+1]);
  
  %switch ib
  %  case 1,
  %    subplot(2,6,[1,2.3]);
  %  case 2,
  %    subplot(2,6,[3.65,4.95]);
  %  case 3,
  %    subplot(2,6,[7,8.3]);
  %  case 4,
  %    subplot(2,6,[9.65,10.95]);
  
  switch icont
    case 1,
      subplot('Position',[0.1 0.66 0.22 0.29]);
    case 2,
      subplot('Position',[0.415 0.66 0.22 0.29])
    case 3,
      subplot('Position',[0.73 0.66 0.22 0.29])
    %case 4,
    %  subplot('Position',[0.1 0.33 0.2 0.23])
    %case 5,
    %  subplot('Position',[0.43 0.33 0.2 0.23])
  end

  axis ([xmin xmax ymin ymax]);
  %axis square;

   % box
  
  ylabel('\epsilon_A');
  xlabel('\epsilon_B');
  hold all;
  
  [Ndir,tt]=size(dirs);

  %Ndir=1
  for idir=1:Ndir
    %file=['../modif/' dirs(idir,1:dirlon(idir)) '_101_' num2str(stgv(i)) '.dat'];
    file=['../modif_B/' dirs(idir,1:dirlon(idir)) '_' bnam '_' num2str(stgv(i)) '_101.dat'];
    disp(file)
    Cfile=importdata(file,' ',1);
    [Ny,Nb]=size(Cfile.data);
    Ny=Ny-1;
    Nb=(Nb-1)/2;
    XX=Cfile.data(2:Ny+1,1);
    YY=Cfile.data(2:Ny+1,ib+1);
    YYinv=Cfile.data(2:Ny+1,Nb-ib+2);
    Zhead=Cfile.data(1,2:Nb+1);
    title (['y=' num2str(Zhead(ib))]);

    
    [xIN,yIN]=polyxpoly(XX,YY,YYinv,XX);
    if length(xIN)==0
      if YYinv(1)>0.8
	iff=find(YY==0); xIN=XX(iff(1)); yIN=0;
      elseif YY(length(YY))>0.8
        xIN=0; iff=find(YYinv==0); yIN=XX(iff(1));
      end
    else
    


    end
    
    length(xIN)
    length(yIN)
    
    Xpol1=XX(1:length(find(XX<xIN)));Ypol1=YY(1:length(find(XX<xIN)));
    Xpol2=flipud(YYinv(1:length(find(XX<yIN)))); Ypol2=flipud(XX(1:length(find(XX<yIN))));
    
    Xpol=cat(1,Xpol1,xIN);Ypol=cat(1,Ypol1,yIN);
    Xpol=cat(1,Xpol,Xpol2);Ypol=cat(1,Ypol,Ypol2);
    
    aa=-(1-Zhead(ib))/(Zhead(ib)); bb=yIN-aa*xIN;
    
    tbb=Ypol-Xpol*aa;
    [tbmax,ibmax]=max(tbb);Xpolmax=Xpol(ibmax); Ypolmax=Ypol(ibmax);
    h=refline([aa tbmax]); set(h,'Color','k');  
    
    Xpol=[0]; Ypol=[0];
    Xpol=cat(1,Xpol,0);Ypol=cat(1,Ypol,tbmax);
    Xpol=cat(1,Xpol,-tbmax/aa);Ypol=cat(1,Ypol,0);
    %plot(Xpol,Ypol,'y');
    h=fill(Xpol,Ypol,'y'); set(h,'EdgeColor','none');   
    
    tbred=inter(2)-inter(1)*aa;
    [xINred,yINred]=polyxpoly(XX,YY,[0 inter(1) 1],[tbred inter(2) aa+tbred]);
    if ( (xINred>xIN) ) 
      xINred=xIN; yINred=yIN;
      tbred=yINred-xINred*aa;
    end
    if ( (inter(1)>=1) )
      xINred=xIN; yINred=yIN;
      tbred=yINred-xINred*aa;
    end    
    if (xINred>0)   
      h=fill([XX(XX<xINred); xINred; (1-tbred)/aa; 0],[YY(XX<xINred); yINred; 1; 1],[1,0.5,0]);  set(h,'EdgeColor','none'); 
    end
  
    
    Xpol=[0]; Ypol=[0]; 
    Xpol=cat(1,Xpol,Xpol1);Ypol=cat(1,Ypol,Ypol1);
    Xpol=cat(1,Xpol,xIN);Ypol=cat(1,Ypol,yIN);
    Xpol=cat(1,Xpol,Xpol2);Ypol=cat(1,Ypol,Ypol2);
    %plot(Xpol,Ypol,'g');
    h=fill(Xpol,Ypol,'g'); set(h,'EdgeColor','none');

 
    h1=plot(XX,YY);
    h2=plot(YYinv,XX);
    set(h1,'Color', colordir(1,:));
    set(h2,'Color', colordir(1,:));
    arx=zeros(Ny,1); ary=zeros(Ny,1)+0.1;
    h=quiver(XX(1:3:length(XX)),YY(1:3:length(YY)),arx(1:3:length(XX)),ary(1:3:length(YY)),0.3);
    set(h,'Color', colordir(1,:), 'ShowArrowHead','off');
    h=quiver(YYinv(1:3:length(YYinv)),XX(1:3:length(XX)),ary(1:3:length(YYinv)),arx(1:3:length(XX)),0.3);
    set(h,'Color', colordir(1,:), 'ShowArrowHead','off');
    %arx=zeros(Ny,1)+0.1;
    %ary=zeros(Ny,1)-0.1*(1-Zhead(ib));
    %h=quiver(XX,YY,arx,ary,0.7);
    %set(h,'Color', 'r');
    %h=quiver(YYinv,XX,-arx,-ary,0.7);
    %set(h,'Color', 'r');
    
    plot([0 1], [0 1], 'k:');

    
    xlim([xmin xmax]);  ylim([ymin ymax]);
    
  
    
  end

  axis ([0.01 0.99 0.01 0.99]);


  
  box
  
  set(gca,'Layer','top')
  
  %rectangle ('Position',[0.01 0.01 0.98 0.98]);
  %set(gca, 'LineWidth',  2)

end


function [xr,yr,ta,tb]=tangc(x,y)
  for i=2:length(x)-2
    xr(i-1)=x(i);
    yr(i-1)=y(i);
    ta(i-1)=(y(i-1)-y(i+1))/(x(i-1)-x(i+1));
    tb(i-1)=yr(i-1)-xr(i-1)*ta(i-1);
  end
end
