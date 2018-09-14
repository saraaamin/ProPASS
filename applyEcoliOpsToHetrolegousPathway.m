% this function is to apply the operatos on two different scenarios:
% Scenario 1:
%   apply the operators associated with the enzymes of the non-native
%   reaction to all the metabolites in E. coli with high concentration
% Scenario 2:
%   apply the operators associated to all the reactions in E.coli on the 
%   main substrates/products of the non-native reactions  
function applyEcoliOpsToHetrolegousPathway()
clear 

% list of global variable that will be used throughout the tool
global reactions
global indexReactions
global compounds
global indexCompounds
global reactionData
global EColiKeggID
global mainModel

% The updated set of data files extracted from KEGG on Jan-27-2016
load reactions_012716
load indexReactions_012716
load compounds_012716.mat
load indexCompounds_012716.mat
load reactionData_012716
load EcoliiJO1366_inGlcOx.mat
load EcoliiJO1366_inGlcOx_KEGGID.mat

% mainModel is the EColi model
originalMainModel = mainModel;

% Keeping track of the Kegg IDs for metabolites in the E.coli model
originalEColiKeggID = EColiKeggID;

% compound is the target metabolite's KEGG ID
compound = 1013;
% compound = 3044;
% compound = 6142;

% biomassPos is biomass reaction number in the model
biomassPos = 7; % Ecoli iJO1366;
biomassMax = 0.9824; %Ecoli iJO1366;

% Files containing the non-native pathways to be studied for different
% target metabolites
load trees_1013
% load trees_3044
% load trees_6142.mat

couplePercent = 0.001;
scenarioFlag = 1;

% this flag pick the set of reactions depending on which scinario is
% desired
if scenarioFlag == 1
    modelRxnIDList = treesUnique{1,1};
else
    modelRxnIDList = mainModel.rxnKEGGIDs;
end

for pathwayIdx = 1:length(treesUnique)
    pathwayIdx
    
    % adding non-native pathway to the metabolic model
    % calculate the flux value for the non-native metabolite through
    % non-native pathway so it would be a point of comparision later on
    [updatedModel, ecoliKeggID, selectedPaths, selectedPathsFlux, wildTypeFluxValue, rxnCompoundList] = CalculatingFluxCell(originalMainModel, originalEColiKeggID, compound,...
        randSelectedPaths(1,pathwayIdx), biomassPos, biomassMax, metsInpathListUnique{pathwayIdx}) ;
    
    % metabolic model with added non-native pathway
    mainModel = updatedModel;
    EColiKeggID = ecoliKeggID;
    
    % retrieve list of enzymes catalyzing the reactions
    enzymesList = getEnzymeIDInEcoli(modelRxnIDList);
    
    % replace EC number of reaction R03544 with 1.1.1.1 maunally for butanol 
    % test case since the one returned by the funtion is too general and not
    % recognized in the operators list
    % pathwayFluxDetails = getSubstrateProductInEcoliRxns(highFluxRxnKeggIDs, enzymesList, EColiKeggID);

    
    % This is the function used to call PROXIMAL and all related functions
    [selectedOperators, allProductsDetails, substrateIDs, prodIDsList, rxnIDList, stepsFBAResults, stepsFBAResultsPos_OneStep, stepsFBAResultsNeg_OneStep,...
        originalRxn_couplingList, allModelPos_OneStep, allModelNeg_OneStep, allModel_Combined, operatorsPerStepList,...
        cids, cidNames, unknownProdIDs, unknownFormulasList] =...
        runProximal_1(enzymesList, biomassPos, couplePercent, metsInpathListUnique{pathwayIdx}, modelRxnIDList, scenarioFlag);
    EColiKeggID_WithOperators{pathwayIdx, 1} = EColiKeggID;
    
    mutantFlux = optimizeCbModel(mainModel, 'max', false, false);

    [equationsList, equationsIDList, formulasList, mainModel] = writeChemicalEqs(allModel_Combined, EColiKeggID, rxnIDList, cids, cidNames);
    balancedMutModel = optimizeCbModel(mainModel, 'max', false, false);    
end

if scenarioFlag == 1
    save mainModel_1013_operators treesUnique selectedPathsFlux EColiKeggID_WithOperators mutantFlux balancedMutModel...
    allProductsDetails substrateIDs prodIDsList rxnIDList stepsFBAResults stepsFBAResultsPos_OneStep stepsFBAResultsNeg_OneStep originalRxn_couplingList...
    allModelPos_OneStep allModelNeg_OneStep allModel_Combined operatorsPerStepList equationsList equationsIDList formulasList selectedOperators...
    cids cidNames unknownProdIDs unknownFormulasList %mainModel_WithOperators_biomassFlux
else
    save mainModel_1013_operators_nonNativeMets treesUnique selectedPathsFlux EColiKeggID_WithOperators mutantFlux balancedMutModel...
    allProductsDetails substrateIDs prodIDsList rxnIDList stepsFBAResults stepsFBAResultsPos_OneStep stepsFBAResultsNeg_OneStep originalRxn_couplingList...
    allModelPos_OneStep allModelNeg_OneStep allModel_Combined operatorsPerStepList equationsList equationsIDList formulasList selectedOperators...
    cids cidNames unknownProdIDs unknownFormulasList %mainModel_WithOperators_biomassFlux
end


