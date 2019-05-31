%Pretreatment of a data set to bild cubic splines input vector

in_v=b2;
clear in_v1 thrs0 t thrs1 n_l
t_simul=40;%%120;%%   total simulated period in seconds
tstep=t_simul/(size(in_v,2));%considering a maximun duration of 40sec

%=======================================================================
i=0;c1=0;
n_l=mean(in_v);%standard noise level in the original signal
clear spk_v0 spk_t0 r
while(i<size(in_v,2)-1)
    i=i+1;        
    if(in_v(i)>n_l)        
        r=1;
        clear a_v
        a_v(r)=in_v(i);        
        while(in_v(i+r)>n_l&&i+r<size(in_v,2))
            r=r+1;
            a_v(r)=in_v(i+r-1);            
        end
        [mx,idx]=max(a_v);
        c1=c1+1;
        spk_v0(c1)=mx;
        spk_t0(c1)=i+idx-1;
        if(r>1)i=i+r;end        
    end
end



%=======================================================================
%this part determines the threshold time_lag assumed periodicity of the
%signal
c1=0;
tOld=0;
for(i=1:size(spk_t0,2))    
        c1=c1+1;
        deltat_V(c1)=spk_t0(i)-tOld;
        tOld=spk_t0(i);
    
end

clear bin_n bin_v h_v
bin_n=200;                          %bins number
bstep=max(deltat_V)/bin_n;          %step between bins
bin_v=bstep/2:bstep:max(deltat_V);  %binned vector
h_v=hist(deltat_V,bin_n)/size(deltat_V,2);         

thrs1=6% 
%=======================================================================
%This part select the participating elements for the construction of the
%splines coefficients
clear deltat spk_v spk_t 
jitter=0.5;
c1=0;
t=1;
%while(in_v(t)==0)t=t+1;end
while(t<=size(in_v,2))    
    if(in_v(t)>0)        
       clear a_v;
       r1=0; r0=t-floor(jitter*thrs1);
       if(r0<=0)r0=t;end
       for r=r0:t+ceil(jitter*thrs1);
            if(r<size(in_v,2))
                r1=r1+1;                              
                a_v(r1)=in_v(r); 
            end
        end
        [mx,idx]=max(a_v);
        c1=c1+1;
        spk_v(c1)=mx;
        spk_t(c1)=tstep*(t-floor(jitter*thrs1)+idx-1);
        t=t-floor(jitter*thrs1)+idx+thrs1-1;
    else
        t=t+1;
    end
        
end

t_v1=tstep:tstep:t_simul;
% % %% ++++++++++++++++++++computing the interpolated values++++++++++++++++++++
splines_V=spline(spk_t,spk_v);
y_val=ppval(splines_V,t_v1);
perc_thr2=.25;
%%=================================================================
%%===========This part computes the waxing/waning histograms=======
 flag_t=2;  dt_u=0;     dt_d=0; i=0; bin_n2=100; bin_n3=100;
 t_old0=0;   t_old1=0;   %%thr2=median(y_val);%%waxing criterion 
clear time_down time_up

if(idx_1==0)thr2=perc_thr2*80;%% 25% of total number of cells
elseif(idx_1==1)thr2=perc_thr2*20;
end
 while(i<size(y_val,2)-1)
     i=i+1;
    switch flag_t        
        case 0                       
            if(y_val(i+1)>thr2&&i<size(y_val,2))%% time down counter                
                dt_d=dt_d+1;
                time_down(dt_d)=t_v1(i)-t_old0;
                flag_t=1;
                t_old1=t_v1(i);
            elseif(y_val(i+1)<=thr2&&i+1==size(y_val,2)&&dt_u>0)%%si llega al final y no cambia
                dt_d=dt_d+1;
                time_down(dt_d)=t_v1(i)-t_old0;    
            elseif(y_val(i+1)<=thr2&&i==size(y_val,2)-1&&dt_u==0)%%fully non synch firing
              dt_d=dt_d+1;
              time_down(dt_d)=t_v1(i)-t_old0;  
              dt_u=dt_u+1;             
              time_up(dt_u)=0; 
            end            
            
                       
            
        case 1  %Waxing phase
            if(y_val(i+1)<=thr2)%%time up condition
                dt_u=dt_u+1;
                time_up(dt_u)=t_v1(i)-t_old1;
                flag_t=0;
                t_old0=t_v1(i);               
            elseif(y_val(i+1)>thr2&&i+1==size(y_val,2)&&dt_d>0)
                dt_u=dt_u+1;
                time_up(dt_u)=t_v1(i)-t_old1;        
            elseif(y_val(i+1)>thr2&&i==size(y_val,2)-1&&dt_d==0)%%fully synch firing
              dt_u=dt_u+1;
              time_up(dt_u)=t_v1(i)-t_old1;  
              dt_d=dt_d+1;             
              time_down(dt_d)=0;              
            end                      
            
             
             
             
        case 2               
            if(y_val(i)>thr2)flag_t=1; t_old1=t_v1(i);
            elseif(y_val(i)<=thr2)flag_t=0; t_old0=t_v1(i); 
            end                
        end
 end


%%==================================================
%%==================================================
mn=0; 
for mm=1:size(spk_t,2)
if(spk_v(1,mm)>0&&spk_t(1,mm)>=0.2)
    mn=mn+1;
    spk_tA(1,mn)=spk_t(1,mm);
    spk_vA(1,mn)=spk_v(1,mm); 
    
end
end
 %%==================================================
 %%==================================================
 
exst=exist('time_up','var'); 
if(exst==0&&flag_t==1)time_up2=log10(t_simul);time_down2=0;elseif(exst==0&&(flag_t>1||flag_t<1))time_up2=log10(t_simul);time_down2=0;; else time_up2=log10(time_up); 
fprintf(1,'the mean value for up is  %f+/-%f \t',mean(time_up),std(time_up))%
end
exst=exist('time_down','var');exst1=exist('time_up2','var');
if(exst1==1 && exst==1)time_down2=log10(time_down);
fprintf(1,'and for down %f+/-%f \n',mean(time_down),std(time_down))%
end
ext2=exist('time_up');
if(ext2==1)
    [hist_HAE,time_bins]=hist(time_up,0:40/1200:40);
else
    time_up=0;
    [hist_HAE,time_bins]=hist(time_up,0:40/1200:40);
end

%%========================================================
%%========================================================
%%plotting bars

bar(t_v1,in_v,'k')
hold on
plot(t_v1,y_val,'r')
set(gca,'FontSize',13);
xlim([10 16]); ylim([0 90]);
set(gca,'FontSize',13,'XTick',[10,11,12,13,14,15,16],'YTick',0:20:80,'ticklength',1.5*get(gca,'ticklength'),'TickDirMode','manual','TickDir','out')
ylabel('Spikes/6ms bin','FontSize',14); xlabel('time(s)','FontSize',14);