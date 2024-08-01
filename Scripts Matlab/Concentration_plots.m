%%% Run gráficos Concentración
file1 = readtable("Ref.txt");
time = 60*file1.Var1-60*20;
Acetref = 16.93*file1.Var2;

file2 = readtable("Mut.txt");
time2 = 60*file2.Var1-60*20;
Acmut = 16.93*file2.Var2;

file3=readtable("Mls.txt");
time3 = 60*file3.Var1-60*20;
Acmls = 16.93*file3.Var2;

file4=readtable("Acs.txt");
time4 = 60*file4.Var1-60*20;
Acs = 16.93*file4.Var2;

%% Concentration plot of 4 experiments
% plot(time,Acetref,"ro",time2,Acmut,"bo",time3,Acmls,"go",time4,Acs,"ko")


 for i=1:45
    if rem(i,9)==0
        subplot(3,3,9);
        plot(tiempo,mets(:,i)/60+20,"r-",tiempo_mut/60+20,mets_mut(:,i),"g-",tiempo_mut2/60+20,mets_mut2(:,i),"b-",tiempo_mut3/60+20,mets_mut3(:,i),"k-")
        grid on
        xlabel("[hr]")
        ylabel('['+nombres(i)+'] mM')
        figure()
    else
        if i==12
            subplot(3,3,rem(i,9))
        plot(tiempo/60+20,mets(:,i),"r-",tiempo_mut/60+20,mets_mut(:,i),"g-",tiempo_mut2/60+20,mets_mut2(:,i),"b-",tiempo_mut3/60+20,mets_mut3(:,i),"k-")
        grid on
        xlabel("[hr]")
        ylabel('['+nombres(i)+'] mM')
        else
           subplot(3,3,rem(i,9))
        plot(tiempo,mets(:,i),"r-",tiempo_mut,mets_mut(:,i),"g-",tiempo_mut2,mets_mut2(:,i),"b-",tiempo_mut3,mets_mut3(:,i),"k-")
        grid on
        xlabel("[hr]")
        ylabel('['+nombres(i)+'] mM')
        end
        
    end
end
subplot(2,1,1)
plot(tiempo/60+20,mets(:,12),"r-",time/60+20,Acetref,"ro",tiempo_mut/60+20,mets_mut(:,12),"g-",time2/60+20,Acmut,"go",tiempo_mut2/60+20,mets_mut2(:,12),"b-",time3/60+20,Acmls,"bo",tiempo_mut3/60+20,mets_mut3(:,12),"k-",time4/60+20,Acs,"ko")
legend("Wild type sim","Wild type","\Delta Cit2 sim","\Delta Cit2 ","\Delta Mls1 sim","\Delta Mls1","\Delta Acs1 sim","\Delta Acs1")
grid on
xlabel("[hr]")
ylabel('[Acetatex] mM')
subplot(2,1,2)
plot(tiempo/60+20,mets(:,1),"r-",tiempo_mut/60+20,mets_mut(:,1),"g-",tiempo_mut2/60+20,mets_mut2(:,1),"b-",tiempo_mut3/60+20,mets_mut3(:,1),"k-")
grid on
xlabel("[hr]")
ylabel('Etoh mM')