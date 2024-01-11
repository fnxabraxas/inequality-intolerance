
function plotlegend()
%  subplot('Position',[0.67 0.23 0.22 0.29])
  subplot('Position',[0.73 0.23 0.22 0.29]);
  axis off;
  axis([0 1 0 1]);
  hold all;
  colordir=[ 0 0.3 1; 0.8 0.03 0.03; 1 1 0; 0 0 0 ];
  transpdir=[1 0.5 0.5 0.5];
%  r=0.25; drh=0.2; ofx=0.09; ofy=0; rr=1-2*r-drh;
%  t = 0:.01:2*pi; xcirc=r*sin(t); ycirc=r*cos(t);
%  fill(xcirc+r,ycirc+0.6,colordir(1,:),'FaceAlpha',transpdir(1),'EdgeColor','none');
%  fill(xcirc+1-r-drh,ycirc+0.6,colordir(2,:),'FaceAlpha',transpdir(2),'EdgeColor','none'); 
%  fill(xcirc+(1-drh)/2,ycirc+0.6-(rr*(3^0.5)/2),colordir(3,:),'FaceAlpha',transpdir(3),'EdgeColor','none');
%  plot(xcirc+r,ycirc+0.6,'k-');
%  plot(xcirc+1-r-drh,ycirc+0.6,'k-'); 
%  plot(xcirc+(1-drh)/2,ycirc+0.6-((1-2*r-drh)*(3^0.5)/2),'k-'); 

  drh=0.2; r=(1-drh)/3; ofx=0.09; ofy=0;
  t = 0:.01:2*pi; xcirc=r*sin(t); ycirc=r*cos(t);
  fill(xcirc+r,ycirc+0.6,colordir(1,:),'FaceAlpha',transpdir(1),'EdgeColor','none');
  fill(xcirc+1-r-drh,ycirc+0.6,colordir(2,:),'FaceAlpha',transpdir(2),'EdgeColor','none'); 
  fill(xcirc+(1-drh)/2,ycirc+0.6-(r*(3^0.5)/2),colordir(3,:),'FaceAlpha',transpdir(3),'EdgeColor','none');
  plot(xcirc+r,ycirc+0.6,'k-');
  plot(xcirc+1-r-drh,ycirc+0.6,'k-'); 
  plot(xcirc+(1-drh)/2,ycirc+0.6-((1-2*r-drh)*(3^0.5)/2),'k-'); 

  
  %rectangle('Position',[0.5-drh-r/2-ofx,0.67-r/2-ofy,r,r],'Curvature',[1,1],'FaceColor',colordir(1,:));
  %rectangle('Position',[0.5+drh-r/2-ofx,0.67-r/2-ofy,r,r],'Curvature',[1,1],'FaceColor',colordir(2,:));
  %rectangle('Position',[0.5-r/2-ofx,0.4-r/2-ofy,r,r],'Curvature',[1,1],'FaceColor',colordir(3,:));
  text(-0.15,0.8,'I|TT','Color','k','FontWeight','bold','Fontsize',12);
  text(0.77,0.8,'T|IT','Color','k','FontWeight','bold','Fontsize',12);
  text(0.34,0.02,'IT|I','Color','k','FontWeight','bold','Fontsize',12);
  
  
%  recx=0.35; drect=0.1; dtrec=0.2; ysep=0.17; ytex=0.93;
%  sizrecx=0.3; sizrecy=0.1; fonts=12;
%  rectangle('Position',[recx,ytex-sizrecy/2,sizrecx,sizrecy],'FaceColor',[0 0.3 1]);
%  rectangle('Position',[recx,ytex-ysep-sizrecy/2,sizrecx,sizrecy],'FaceColor',[246 142 142]/255);
%  rectangle('Position',[recx,ytex-2*ysep-sizrecy/2,sizrecx,sizrecy],'FaceColor',[250 250 0]/255);
%  rectangle('Position',[recx,ytex-3*ysep-sizrecy/2,sizrecx,sizrecy],'FaceColor',[100 170 60]/255);
%  rectangle('Position',[recx,ytex-4*ysep-sizrecy/2,sizrecx,sizrecy],'FaceColor',[110 30 140]/255);
%  rectangle('Position',[recx,ytex-5*ysep-sizrecy/2,sizrecx,sizrecy],'FaceColor',[250 180 70]/255);
%  text(recx+sizrecx+drect,ytex,'I|TT','Color','k','FontWeight','bold','Fontsize',fonts);
%  text(recx+sizrecx+drect,ytex-ysep,'T|IT','Color','k','FontWeight','bold','Fontsize',fonts);
%  text(recx+sizrecx+drect,ytex-2*ysep,'TT|I','Color','k','FontWeight','bold','Fontsize',fonts);
%  text(recx+sizrecx+drect,ytex-3*ysep,'I|TT + TT|I','Color','k','FontWeight','bold','Fontsize',fonts);
%  text(recx+sizrecx+drect,ytex-4*ysep,'I|TT & T|IT','Color','k','FontWeight','bold','Fontsize',fonts);
%  text(recx+sizrecx+drect,ytex-5*ysep,'T|IT & TT|I','Color','k','FontWeight','bold','Fontsize',fonts);
end