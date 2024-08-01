%%% Run gráficos flujos
file1 = readtable("Ref.txt");
time = 60*file1.Var1-60*20;
Acetref = 16.93*file1.Var2;
file2 = readtable("Mut.txt");
time2 = 60*file2.Var1-60*20;
Acmut = 16.93*file2.Var2;

 figure()
 for i=1:36
     if rem(i,9)==0
        subplot(3,3,9)
        plot(tiempo/60+20,flujos(:,i),"r-",tiempo_mut/60+20,flujos_mut(:,i),"g-",tiempo_mut2/60+20,flujos_mut2(:,i),"b-",tiempo_mut3/60+20,flujos_mut3(:,i),"k-")
        grid on
         xlabel("[hr]")
        ylabel(''+NOMBRE_flujos(i)+" [(mM/min)]")
        figure()
     else
        subplot(3,3,rem(i,9))
        plot(tiempo/60+20,flujos(:,i),"r-",tiempo_mut/60+20,flujos_mut(:,i),"g-",tiempo_mut2/60+20,flujos_mut2(:,i),"b-",tiempo_mut3/60+20,flujos_mut3(:,i),"k-")
        if i==4
        plot(tiempo/60+20,flujos(:,i),"r-",tiempo_mut/60+20,flujos_mut(:,i),"g-",tiempo_mut2/60+20,flujos_mut2(:,i),"b-",tiempo_mut3/60+20,flujos_mut3(:,i),"k-")
        grid on
         xlabel("[hr]")
        ylabel(''+NOMBRE_flujos(i)+' [(mM/min)]')
        ylim([-0.01 0.04])
        else
            plot(tiempo/60+20,flujos(:,i),"r-",tiempo_mut/60+20,flujos_mut(:,i),"g-",tiempo_mut2/60+20,flujos_mut2(:,i),"b-",tiempo_mut3/60+20,flujos_mut3(:,i),"k-")
            grid on
         xlabel("[hr]")
        ylabel(''+NOMBRE_flujos(i)+' [(mM/min)')

        end
           
     end
 end
 plot(tiempo/60+20,flujos(:,4)/(P.Vtot-P.Vcels),"r-",time/60+20,gradient(Acetref,time),"ro",tiempo_mut/60+20,flujos_mut(:,4)/(P.Vtot-P.Vcels),"g-",time2/60+20,gradient(Acmut,time2),"go",   tiempo_mut2/60+20,flujos_mut2(:,4)/(P.Vtot-P.Vcels),"b-",time3/60+20,gradient(Acmls,time3),"bo",   tiempo_mut3/60+20,flujos_mut3(:,4)/(P.Vtot-P.Vcels),"k-",time4/60+20,gradient(Acs,time4),"ko")
 grid on
 legend("Wild type sim","Wild type","\Delta Cit2 sim","\Delta Cit2 ","\Delta Mls1 sim","\Delta Mls1","\Delta Acs1 sim","\Delta Acs1")
 xlabel("[hr]")
 ylabel(''+NOMBRE_flujos(4)+' [(mM/min)]')