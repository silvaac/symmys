function [A,b,g]=CondProbViews(View,X)
% input views
% statement: View(k).Who (e.g. [1 3])= View(k).Equal (e.g. {[2 3] [1 3 5]})
% optional conditional statement: View(k).Cond_Who (e.g. [2])= View(k).Cond_Equal (e.g. {[1]})
% amount of stress is quantified as Prob(statement) <= View(k).v if View(k).sgn = 1;
%                                   Prob(statement) >= View(k).v if View(k).sgn = -1;
% confidence in stress is quantified in View(k).c in (0,1)

A=[];
b=[];
g=[];
for k=1:length(View)
    
    I_mrg=(X(:,1)<Inf);
    for s=1:length(View(k).Who)
        Who=View(k).Who(s);
        Or_Targets=View(k).Equal{s};
        I_mrg_or=(X(:,Who)>Inf);
        for i=1:length(Or_Targets)
            I_mrg_or = I_mrg_or | (X(:,Who)==Or_Targets(i));
        end
        I_mrg = I_mrg & I_mrg_or;
    end

    I_cnd=(X(:,1)<Inf);
    for s=1:length(View(k).Cond_Who)
        Who=View(k).Cond_Who(s);
        Or_Targets=View(k).Cond_Equal{s};
        I_cnd_or=(X(:,Who)>Inf);
        for i=1:length(Or_Targets)
            I_cnd_or = I_cnd_or | X(:,Who)==Or_Targets(i);
        end
        I_cnd = I_cnd & I_cnd_or;
    end
    
    I_jnt=I_mrg & I_cnd;
    
    if ~isempty(View(k).Cond_Who)
        New_A=View(k).sgn*(I_jnt-View(k).v*I_cnd)';
        New_b=0;
    else
        New_A=View(k).sgn*I_mrg';
        New_b=View(k).sgn *View(k).v;
    end
    
    A = [A
        New_A];  % constraint for the conditional expectations...
    b = [b
        New_b];
    g = [g
        -log(1-View(k).c)];
end
