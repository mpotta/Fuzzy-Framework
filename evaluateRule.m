Xin = csvread('LUNG_mRmR_DAT.csv');
Xout = csvread('LUNG_output.csv');

patterns = size(Xin, 1);
attributes = size(Xin, 2);
candidateNo = 10;

for i = 1:size(Xin,2)
    minValue = min(Xin(:,i));
    maxValue = max(Xin(:,i));
    Xin(:,i) = (Xin(:,i) - minValue)/(maxValue - minValue);
end;

allPermutations = generatePermutations();
compatibility = zeros(size(allPermutations,1),1);
weight = zeros(size(allPermutations,1),1);
myC1Rules = [];
myC2Rules = [];

mftype = 'gaussmf';
mfparams = [[0.5 0]; [0.5 1]; [0.2123 0.5]; [0.2123 0.0]; [0.2123 1.0]];

for k=1:size(allPermutations,1)

    for i=1:patterns        
        compatibilty = ones(patterns,1);
        term = 1;
        
        for j=1:attributes
            term = evalmf(Xin(i,j), mfparams(allPermutations(k,j),:), mftype);
            compatibility(i,1) = compatibility(i,1)*term;
        end;
    end;
    
    for i=1:patterns        
        confidenceC1 = 0;
        confidenceC2 = 0;
        totalCompatibility = 0;
        
        if(Xout(i,1) == 0)
            confidenceC1 = confidenceC1 + compatibility(i,1);
        else
            confidenceC2 = confidenceC2 + compatibility(i,1);
        end;
        
        totalCompatibility = totalCompatibility + compatibility(i,1);
    end;
    
    confidenceC1 = confidenceC1/totalCompatibility;
    confidenceC2 = confidenceC2/totalCompatibility;
    maxConfidence = max(confidenceC1, confidenceC2);
    minConfidence = min(confidenceC1, confidenceC2);
    weight(k,1) = maxConfidence - minConfidence;
    
    if(maxConfidence == confidenceC1)
        myC1Rules(k,:) = [allPermutations(k,:), 0];
    else
        myC2Rules(k,:) = [allPermutations(k,:), 1];
    end;   
end;

sortIndex = size(myC2Rules,2);
myCRules = [myC1Rules; myC2Rules];
myCRules = sortrows(myCRules, -sortIndex);

candidateRules = myCRules(1:candidateNo,:);

% Initialize

Npop = 50;                      % Population Size
Q = size(candidateRules,1);                % Chromosome Length
nMax = 100;                     % Maximum Allowable Generations
ruleProbability = 0.4;          % Probability of Rule Selection

myPopulation = zeros(Npop,Q);
interimPopulation = zeros(Npop,Q);
newPopulation = zeros(Npop,Q);

myFitness = zeros(1,Npop);
interimFitness = zeros(1,Npop);

pc = 0.7;                       % Crossover Probability
pA = 0.05;
w1 = 1.0;
w2 = 0;
w3 = 0;
pA10 = (w2/(w1+w3))*pA;         % Mutation Probability 10
pA01 = ((w1+w3)/w2)*pA;         % Mutation Probability 01

% Initialize Population
for i = 1:Npop
    myPopulation(i,:) = randsrc(1,Q,[0 1; (1-ruleProbability) ruleProbability]);
end;

% Generate Next Generation
for i = 1:nMax 
    % Genetic Operations: Crossover & Mutation    
    for j = 1:Npop
        
        % Two Point Crossover
        parent1 = j;
        parent2 = randperm(Npop,1);
        while isequal(parent1, parent2)
            parent2 = randperm(Npop,1);
        end;
        parent1 = myPopulation(parent1,:);
        parent2 = myPopulation(parent2,:);
        
        if(rand < pc)
                       
            temp1 = randperm(Q,1);
            temp2 = randperm(Q,1);
            while (abs(temp1 - temp2) < 2)
                temp2 = randperm(Q,1);
            end;
            crossOverPoint1 = min(temp1,temp2);
            crossOverPoint2 = max(temp1,temp2);
            
            child = parent1;                        
            child(1,crossOverPoint1:crossOverPoint2-1) = parent2(1,crossOverPoint1:crossOverPoint2-1);
            
        else
            child = parent1;
           
        end;
       
        % Mutation
        for a = 1:Q
            if(child(1,a) == 0)
                if(rand < pA01)
                    child(1,a) = 1;
                end;
            else
                if(rand < pA10)
                    child(1,a) = 0;
                end;
            end;
        end;
        
        interimPopulation(j,:) = child(1,:);
        interimFitness(1,j) = getFitness(child(1,:), candidateRules, w1, w2, w3);
    end;
    
    % Determine New Population
    myPopIndex = [myFitness; 1:Npop];
    myPopIndex = sortrows(myPopIndex')';
    
    interimPopIndex = [interimFitness; 1:Npop];
    interimPopIndex = sortrows(interimPopIndex')';
    
    for k1 = 1:(Npop/2)
        newPopulation(k1,:) = myPopulation(myPopIndex(2,k1),:);
    end;
    
    for k2 = 1:(Npop/2)
        newPopulation((Npop/2)+k2,:) = interimPopulation(interimPopIndex(2,k2),:);
    end;
    
    myPopulation = newPopulation;
    
end;

x=4;

[~, index] = max(myFitness);
solution = myPopulation(index,:);

patterns = size(Xin, 1);
q = size(solution, 2);

accuracy = 0;

a = newfis('fuzzyModel');
for i=1:attributes
    a = addvar(a,'input',strcat('gene', num2str(i)),[0 1]); 
end;
a = addvar(a,'output','cancer', [0 1]);

for i=1:attributes
    a = addmf(a,'input',i,'S1','gaussmf', [0.5 0]);
    a = addmf(a,'input',i,'S2','gaussmf', [0.5 1]);
    a = addmf(a,'input',i,'M3','gaussmf', [0.2123 0.5]);
    a = addmf(a,'input',i,'S3','gaussmf', [0.2123 0.0]);
    a = addmf(a,'input',i,'L3','gaussmf', [0.2123 1.0]);
end;

a = addmf(a,'output',1,'C1','trimf', [0 0.25 0.5]);
a = addmf(a,'output',1,'C2','trimf', [0.5 0.75 1]);

a.rule = [];
ruleList = [];
for i=1:Q
    if(solution(1,i) == 1)
        r = [candidateRules(i,:), 1, 1];
        ruleList = [ruleList; r];
    end;
end;
a = addrule(a,ruleList);

for i=1:patterns;
    output = evalfis(Xin(i,:),a);
    if(output < 0.5) 
        output = 0; 
    else output = 1; 
    end;
    if(output == Xout(i,1)) 
        accuracy = accuracy + 1; 
    end;
end;

accuracy = accuracy/patterns;



         
    
    
  
    
    














