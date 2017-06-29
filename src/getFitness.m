function [fitness] = getFitness(myChromosome, candidateRules, w1, w2, w3)

NCCP = 0;
NOR = sum(myChromosome);
NOA = 5;

Xin = csvread('COLON_mRmR_DAT.csv');
Xout = csvread('COLON_output.csv');

patterns = size(Xin, 1);
attributes = size(Xin, 2);
ruleNo = size(candidateRules, 1);
ruleLength = size(candidateRules, 2);

for i = 1:size(Xin,2)
    minValue = min(Xin(:,i));
    maxValue = max(Xin(:,i));
    Xin(:,i) = (Xin(:,i) - minValue)/(maxValue - minValue);
end;

afis = newfis('fuzzyModel');
for i=1:attributes
    afis = addvar(afis,'input',strcat('gene', num2str(i)),[0 1]); 
end;
afis = addvar(afis,'output','cancer', [0 1]);

for i=1:attributes
    afis = addmf(afis,'input',i,'S1','gaussmf', [0.5 0]);
    afis = addmf(afis,'input',i,'S2','gaussmf', [0.5 1]);
    afis = addmf(afis,'input',i,'M3','gaussmf', [0.2123 0.5]);
    afis = addmf(afis,'input',i,'S3','gaussmf', [0.2123 0.0]);
    afis = addmf(afis,'input',i,'L3','gaussmf', [0.2123 1.0]);
end;

afis = addmf(afis,'output',1,'C1','trimf', [0 0.25 0.5]);
afis = addmf(afis,'output',1,'C2','trimf', [0.5 0.75 1]);

afis.rule = [];
ruleList = [];
for i=1:ruleNo
    if(myChromosome(1,i) == 1)
        r = [candidateRules(i,:), 1, 1];
        ruleList = [ruleList; r];
    end;
end;
afis = addrule(afis,ruleList);

for i=1:patterns;
    output = evalfis(Xin(i,:),afis);
    if(output < 0.5) 
        output = 0; 
    else output = 1; 
    end;
    if(output == Xout(i,1)) 
        NCCP = NCCP + 1; 
    end;
end;


fitness = (w1*NCCP) - (w2*NOR) - (w3*NOA);

return;