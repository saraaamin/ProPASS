% running proximal
% input:
%   enzymesList: a list of EC numbers for which we want to run proximal
%   biomassPos: the position of biomass in the metabolic model
%   couplePercent, metsInpath, rxnKEGGIDList, scenarioFlag
% output:
%   mutantFluxValue: the result of running FBA after augmenting the model
%   with the operators and related reactions
%   selectedOperators: a list of selected operators for the enzymes
%   allProductsDetails: the number of products generated for each
%   metabolite in the E.coli model
%   substrateIDs: a list of substrate IDs taking place in the model
%   augmentation
%   prodIDsList: 
function [selectedOperators, allProductsDetails, substrateIDs, prodIDsList, rxnIDList, ...
    stepsFBAResults, stepsFBAResultsPos_OneStep, stepsFBAResultsNeg_OneStep, originalRxn_couplingList,...
    allModelPos_OneStep, allModelNeg_OneStep, allModel_Combined,...
    operatorsPerStepList, cids, cidNames, unknownProdIDs, unknownFormulasList] = ...
    runProximal_1(enzymesList, biomassPos, couplePercent, metsInpath, rxnKEGGIDList, scenarioFlag)

global EColiKeggID
global mainModel
global operators
global operatorsMainAtom
global operatorsMainAtomPos

% load the Operators Neda generated from KEGG
% load KEGGOperators_2016_08_14.mat
% load KEGGOperators_2016_08_14_update.mat
load operators.mat

[~, rxnsNum] = size(enzymesList);
% for each pathway, get the list of enzymes catalyzing the reactions
selectedOperators = struct();
for rxnIdx = 1:rxnsNum
    enzymesOfCurrentRxn = enzymesList(1,rxnIdx);
    enzymesNum = length(enzymesOfCurrentRxn{1,1});
   
     % each enzyme in the enzymesList, get the list of operators associated
     % to it from the allOperators dataset
     if rxnKEGGIDList(rxnIdx) == 1213
         display('here');
     end
    for enzymeIdx = 1:enzymesNum
        selectedOperators = getEnzymeOperators(enzymesOfCurrentRxn{1,1}{1,enzymeIdx}, rxnKEGGIDList(rxnIdx), selectedOperators);
    end
end

selectedOperators = getUniqueOperators(selectedOperators);

[selectedOperators, operatorsMainAtom, operatorsMainAtomPos] = sortRs(selectedOperators);
save SelectedOperators selectedOperators operatorsMainAtom operatorsMainAtomPos

%     [operatorsPathwayDetails, newOperators] = applyOperatorsToPathway('SelectedOperators.mat', pathwayFluxDetails);

[allProductsDetails, substrateIDs, prodIDsList, rxnIDList, stepsFBAResults, stepsFBAResultsPos_OneStep, stepsFBAResultsNeg_OneStep, ...
    originalRxn_couplingList, allModelPos_OneStep, allModelNeg_OneStep, allModel_Combined, operatorsPerStepList, selectedOperators,...
    cids, cidNames, unknownProdIDs, unknownFormulasList] = ...
    applyOperatorsToModel('SelectedOperators.mat', biomassPos, couplePercent, metsInpath, scenarioFlag);
% save SelectedOperators selectedOperators operatorsMainAtom operatorsMainAtomPos

end


% This is a function that returns the operators associated with a specific
% enzyme.
% input:
%   enzyme: EC number of the quired enzyme
%   selectedOperators: a list of current selected operators in case the
%   function is used to find the operators for more than one enzyme
%   catalyzing the same reaction
% output:
%   selectedOperators: a cummulitave list of all operators to be considered
function selectedOperators = getEnzymeOperators(enzymeList, rxnID, selectedOperators)

global operators
operatorsNum = length(operators);

% if a reaction has more than one enzyme associated to it, operators are
% appended in one list, and below is to keep track of the number of items
% in the operators list so far
if isempty(fieldnames(selectedOperators))
    selectedOperatorsNum = 0;
else
    selectedOperatorsNum = length(selectedOperators);
end

% get the operators of the EC number passed to the function from the
% operators dataset, and operator has to include the enzyme and the
% reaction number to be selected. Some operators don't include all the
% reactions related to them in Kegg
for enzymeIdx = 1:length(enzymeList)
    for operatorIdx = 1:operatorsNum
        enzymeMatch = find(strcmp(enzymeList{enzymeIdx}, operators(operatorIdx).Enzyme));
        rxnMatch = find(operators(operatorIdx).Reaction == rxnID);
        
        if ~isempty(enzymeMatch) && ~isempty(rxnMatch)
            selectedOperatorsNum = selectedOperatorsNum+1;
            selectedOperators(selectedOperatorsNum).ID = operators(operatorIdx).ID;
            selectedOperators(selectedOperatorsNum).Reaction = operators(operatorIdx).Reaction;
            selectedOperators(selectedOperatorsNum).Enzyme = operators(operatorIdx).Enzyme;
            selectedOperators(selectedOperatorsNum).Reactant = operators(operatorIdx).Reactant;
            selectedOperators(selectedOperatorsNum).Product = operators(operatorIdx).Product;
            selectedOperators(selectedOperatorsNum).KCF = operators(operatorIdx).KCF;
            
            selectedOperators(selectedOperatorsNum).subKEGGID = operators(operatorIdx).SubstrateID;
            selectedOperators(selectedOperatorsNum).prodKEGGID = operators(operatorIdx).ProductID;
            selectedOperators(selectedOperatorsNum).CoupledRxn = rxnID;   
            
%             selectedOperators(selectedOperatorsNum).CoupledRxn = pathwayDetails.rxn;
%             selectedOperators(selectedOperatorsNum).susbtrateSideCofactors = pathwayDetails.substrateSideCofactors;
%             selectedOperators(selectedOperatorsNum).productSideCofactors = pathwayDetails.productSideCofactors;
        end
    end
end
end

% This is a function that applies a set of operators to each compound in 
% the model provided to generate products to be augmented to the model
% input:
%   operatorsFileName: the name of the operators mat file
%   biomassPos: the position of biomass in the model
%   used to produce a target molecule, and the enzyme and flux value for
%   each reaction
% output:
%   allProductsDetails: an array containing the number of products
%   resulting from applying the operators to each compound in the model
%   substrateIDs: a list of KEGG IDs for compounds that generated products
%   when operators are applied to them
%   prodIDsList: a list of KEGG IDs for the products produced by an
%   operator, they correspond to the sustrateIDs
%   rxnIDList: a list of reactions that are updated or added to the model
%   when a new product is identified
%   stepsFBAResults: a structure containing the FBA results of maximizing 
%   the flux of the a target when when a new product is added to the model
%   originalRxn_couplingList: a list of reactions IDs (from the synthesis 
%   pathway) to which each step agumented to the model is coupled with
function [allProductsDetails, substrateIDs, prodIDsList, rxnIDList, stepsFBAResults, stepsFBAResultsPos_OneStep, stepsFBAResultsNeg_OneStep,...
    originalRxn_couplingList, allModelPos_OneStep, allModelNeg_OneStep, allModel_Combined, operatorsPerStepList, selectedOperators,...
    cids, cidNames, unknownProdIDs, unknownFormulasList]  = ...
    applyOperatorsToModel(operatorsFileName, biomassPos, couplePercent, metsInpath, scenarioFlag)

load cofactors.mat

global EColiKeggID
global mainModel
global substrateIDs
global prodIDsList

originalModel = mainModel;
originalEColiKeggID = EColiKeggID;

rxnIDList = [];
allProductsDetails = [];
productsFolder = '.\productsMol\';
savedProductsFolder = '.\savedProductsMol\';
delete([savedProductsFolder, '\*.mol']);

uniqueEColiKeggID = unique(EColiKeggID);
cids = [];
cidNames = [];
cidCount = 0;
unknownProdIDs = [];
unknownFormulasList = [];
substrateIDs = [];
prodIDsList = [];
stepsFBAResults = [];
stepsFBAResultsPos_OneStep = [];
stepsFBAResultsNeg_OneStep = [];
% stepsFBAResults_biomass = [];
allModelPos_OneStep = [];
allModelNeg_OneStep = [];
allModel_Combined=[];
originalRxn_couplingList = [];
operatorsPerStepList = [];
newProdID = -101.01;

load(operatorsFileName)
load ecoli_concentration_data.mat

if scenarioFlag == 1
   compoundsToApplyOpsOn = uniqueEColiKeggID;
else
   compoundsToApplyOpsOn = metsInpath;  
end
% for each compound in the model, apply the operators to generate products
% associated to the operators

for compoundIdx = 1:length(compoundsToApplyOpsOn)
    compound = compoundsToApplyOpsOn(compoundIdx);

    % this condition is to check if the compound has concentration less
    % than 10^-6 then it will be ignored
    if scenarioFlag == 1 && isempty(find(filteredMets_KEGGIDs == compound))
        continue;
    end
    
%     compound IDs that don't have matches in KEGG or compounds consist of
%     one atoms with no bounds. Those compounds are skipped since no
%     products can be generated from them
%     the last set of IDs are carrier proteins that should be excluded
%     since they contain S and R groups which can't be balanced as they are
%     not detailed.
    excludedIDs = [-1; 0; 1; 10; 23; 34; 38; 70; 76; 80; 84; 87; 175; 238; 282; 283; 291; 305; 698; 703; 787; 824; 1330; 1342; 1413; 1528; 1635;...
        1636; 1637; 1638; 1639; 1640; 1641; 1642; 1643; 1644; 1645; 1646; 1647; 1648; 1649; 1650; 1651; 1652; 1653; 1834; 2386; 2745; 2869;...
        5737; 6710; 7292; 14818; 14819; 15233; 19610;...
        173; 4180; 4619; 4620; 4688; 5223; 5274];
    
    
    if ~isempty(find(excludedIDs == compound))
        continue;
    end
    
    if ~isempty(find(gluCoFactKEGGID == compound))
        continue;
    end
    
    
%     get the string of a compound ID to look up the KCF file on KEGG
    digitsNum = floor(log10 (compound))+1;
    f = '';
    if digitsNum~=5
        f = '0';
        for i = 1:5-1-digitsNum
            f = strcat(f,'0');
        end
    end
    
    compoundstr = strcat('C',f,num2str(compound));
    % generate a list of products when applying the operators to the
    % current compound 
    [inputList, inputStructure, inputListNumbering, products, operatorIdx] = GenerateProducts(compoundstr,operatorsFileName);
    allProductsDetails(compoundIdx) = length(products);
    
    
    
    if ~isempty(products)
        prodEnzymesList = GenerateMolFiles(inputStructure, inputListNumbering, products);
        
        productsFolderFiles = dir(productsFolder);
        
        for prodIdx = 3:length(productsFolderFiles)
            prodIdx    
            prodFormula = [];
            currentProdFile = productsFolderFiles(prodIdx).name;
            fileIdx = str2double(productsFolderFiles(prodIdx).name(9:end-4));
            prodKCFFile = generate_kcfFile([productsFolder , currentProdFile]);
            prodCmpID = compareProdKCFToKeggKCF(prodKCFFile);
            
            if prodCmpID == compound
                continue;
            end
                
            % Since omega-Carboxyacyl-CoA is a general metabolite, it's not
            % considered as a viabal product
            
            % operatorsPathwayDetailsIdx = getReactionCofactors(selectedOperators, operatorIdx(prodIdx-2), prodEnzymesList{prodIdx-2,1}, pathwayFluxDetails);
            % operatorsPathwayDetailsIdx = operatorIdx(prodIdx-2);
            operatorsPathwayDetailsIdx = operatorIdx(fileIdx);
            
            if ~isempty(prodCmpID)
                if length(prodCmpID) == 1
                    switchCase = 1;
                else
                    switchCase = 2;
                end
            else
                % retrieve details using OpenBabel toolbox
                switchCase = 1;
                productOBDetails = py.OpenBabelFilesConversion.main_function(currentProdFile);
                
                if ~isempty(productOBDetails)% && py.len(productOBDetails{1}) > 
                    prodOBStr = char(py.str(productOBDetails));
                    separatorIdx = find(prodOBStr == ',');
                    if ~isempty(separatorIdx)
                        pubChemID = str2double(prodOBStr(1:separatorIdx-1));
                        if ~isnan(pubChemID)
                            prodCmpID = -1*pubChemID;
                            pubChemName = prodOBStr(separatorIdx+1:end);
                            cidCount = cidCount + 1;

                            if isempty(find(prodCmpID == cids))
                                cids(cidCount,1) = prodCmpID;
                                cidNames{cidCount,1} = pubChemName;
                            end
                        else
%                             continue;
                            prodCmpID = newProdID;
                            prodFormula = py.OpenBabelFilesConversion.convert_mol_to_formula(currentProdFile);
                            prodFormula = char(py.str(prodFormula));
                            unknownProdIDs(end+1,1) = prodCmpID;
                            unknownFormulasList{end+1, 1} = prodFormula;
                            newProdID = newProdID - 1;
                        end
                    else
%                         continue;
                        prodCmpID = newProdID;
                        prodFormula = py.OpenBabelFilesConversion.convert_mol_to_formula(currentProdFile);
                        prodFormula = char(py.str(prodFormula));
                        unknownProdIDs(end+1,1) = prodCmpID;
                        unknownFormulasList{end+1, 1} = prodFormula;
                        newProdID = newProdID - 1;
                    end
                   
                else
%                     continue;
                    prodCmpID = newProdID;
                    prodFormula = py.OpenBabelFilesConversion.convert_mol_to_formula(currentProdFile);
                    prodFormula = char(py.str(prodFormula));
                    unknownProdIDs(end+1,1) = prodCmpID;
                    unknownFormulasList{end+1, 1} = prodFormula;
                    newProdID = newProdID - 1;

                end
                
            end  
            
            
            if ~isempty(find(gluCoFactKEGGID == prodCmpID(1)))
                continue;
            end
            
            if isempty(operatorsPathwayDetailsIdx)
                continue;
            end
            
            if switchCase == 1
                index = find(originalEColiKeggID == prodCmpID);
                [updatedRxnIDList, updatedModel, modelAugmentFlag] = updateModel_OneStep(originalModel, originalEColiKeggID, index, compound, prodCmpID, selectedOperators(operatorsPathwayDetailsIdx));

                
                if modelAugmentFlag
                    % couple the reaction so the flux would be in the same
                    % direction as the original reaction flux
%                     [updatedModel_PosRxn, updatedEColiKeggID, mainCouplingRxn_PosRxn] = coupleRxnFlux(updatedModel, originalEColiKeggID, prodEnzymesList{prodIdx-2,1},...
%                         pathwayFluxDetails, updatedRxnIDList(end), couplePercent, 1); 
%                     [updatedModel_PosRxn, updatedEColiKeggID, mainCouplingRxn_PosRxn] = restrictFluxOfAddedReaction(compound, updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, 1, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1});
                    [updatedModel_PosRxn, updatedEColiKeggID, mainCouplingRxn_PosRxn, posBoundVal] = setBoundsForAddedRxn(updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, 1, prodEnzymesList{prodIdx-2,1}, selectedOperators(operatorsPathwayDetailsIdx));
                    mutantPosFluxValue = optimizeCbModel(updatedModel_PosRxn, 'max', false, false);

                    % couple the reaction so the flux would be in the
                    % opposite direction as the original reaction flux
%                     [updatedModel_NegRxn, updatedEColiKeggID, mainCouplingRxn_NegRxn] = coupleRxnFlux(updatedModel, originalEColiKeggID, prodEnzymesList{prodIdx-2,1},...
%                         pathwayFluxDetails, updatedRxnIDList(end), couplePercent, -1);

%                     [updatedModel_NegRxn, updatedEColiKeggID, mainCouplingRxn_NegRxn] = restrictFluxOfAddedReaction(compound, updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, -1, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1});
                    [updatedModel_NegRxn, updatedEColiKeggID, mainCouplingRxn_NegRxn, negBoundVal] = setBoundsForAddedRxn(updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, -1, prodEnzymesList{prodIdx-2,1}, selectedOperators(operatorsPathwayDetailsIdx));
                    mutantNegFluxValue = optimizeCbModel(updatedModel_NegRxn, 'max', false, false);


                    % setting default for copuling to be in the same direction
                    coupleDirection = 1;
                    mutantFluxValue = mutantPosFluxValue; 
                    boundVal = posBoundVal;
                    
                    % check if the default coupling direction will change                
                    if mutantPosFluxValue.f < mutantNegFluxValue.f && mutantPosFluxValue.f > 0.0001
                        coupleDirection = 1;
                        mutantFluxValue = mutantPosFluxValue; 
                        boundVal = posBoundVal;
                    elseif mutantNegFluxValue.f < mutantPosFluxValue.f && mutantNegFluxValue.f > 0.0001
                        coupleDirection = -1;
                        mutantFluxValue = mutantNegFluxValue;
                        boundVal = negBoundVal;
                    end
                    
                    modelAugmentFlag = false;

                    if mutantFluxValue.f > 0.0001
                        index = find(EColiKeggID == prodCmpID);
                        [rxnIDList, modelAugmentFlag] = updateMainModel(index, compound, prodCmpID, rxnIDList, selectedOperators(operatorsPathwayDetailsIdx));
                    end

                    if modelAugmentFlag
                        stepsFBAResultsPos_OneStep = [stepsFBAResultsPos_OneStep; mutantPosFluxValue];
                        allModelPos_OneStep = [allModelPos_OneStep; updatedModel_PosRxn];

                        stepsFBAResultsNeg_OneStep = [stepsFBAResultsNeg_OneStep; mutantNegFluxValue];
                        allModelNeg_OneStep = [allModelNeg_OneStep; updatedModel_NegRxn];

%                         [mainModel, EColiKeggID, mainCouplingRx] = restrictFluxOfAddedReaction(compound, mainModel, EColiKeggID, rxnIDList(end), couplePercent, coupleDirection, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1});
                        [mainModel, EColiKeggID, mainCouplingRx] = setBoundsForAddedRxn_FinalModel(mainModel, EColiKeggID, rxnIDList(end), coupleDirection, prodEnzymesList{prodIdx-2,1}, boundVal, selectedOperators(operatorsPathwayDetailsIdx));
                    
%                         [mainModel, EColiKeggID, mainCouplingRxn] = coupleRxnFlux(mainModel, EColiKeggID,...
%                             prodEnzymesList{prodIdx-2,1}, pathwayFluxDetails, rxnIDList(end), couplePercent, coupleDirection); 
    %                     mutantFluxValue = optimizeCbModel(mainModel, 'max', false, false);
    %                     stepsFBAResults = [stepsFBAResults; mutantFluxValue];
    %                     originalRxn_couplingList = [originalRxn_couplingList; mainCouplingRxn];

                        selectedOperators(operatorsPathwayDetailsIdx).moleculePosition = products(fileIdx).ObtainedFromInput;
                        selectedOperators(operatorsPathwayDetailsIdx).coupleDirection = coupleDirection;
                        selectedOperators(operatorsPathwayDetailsIdx).OperatorIdx = operatorsPathwayDetailsIdx;

                        operatorsPerStepList = [operatorsPerStepList; selectedOperators(operatorsPathwayDetailsIdx)];
                        allModel_Combined = [allModel_Combined; mainModel];

                        % saving the mol file in a different directory for
                        % later consideration
                        newNameProdFile = [currentProdFile(1:end-4), '-', num2str(prodCmpID), '.mol'];
                        sourceFolder = [productsFolder, currentProdFile];
                        destinationFolder = [savedProductsFolder, newNameProdFile];
                        movefile(sourceFolder, destinationFolder);
                        
                    end

                end
                
                
            elseif switchCase == 2
                currentCouplePercent = couplePercent/length(prodCmpID);
                
                for kfcProdIdx = 1:length(prodCmpID)
                    if ~isempty(find(gluCoFactKEGGID == prodCmpID(kfcProdIdx),1))
                        continue;
                    end

                    index = find(originalEColiKeggID == prodCmpID(kfcProdIdx));
                    [updatedRxnIDList, updatedModel, modelAugmentFlag] = updateModel_OneStep(originalModel, originalEColiKeggID, index, compound, prodCmpID(kfcProdIdx), selectedOperators(operatorsPathwayDetailsIdx));
                   
                    if modelAugmentFlag
                        % couple the reaction so the flux would be in the same
                        % direction as the original reaction flux
%                         [updatedModel_PosRxn, updatedEColiKeggID, mainCouplingRxn_PosRxn] = coupleRxnFlux(updatedModel, originalEColiKeggID, prodEnzymesList{prodIdx-2,1},...
%                                             pathwayFluxDetails, updatedRxnIDList(end), couplePercent, 1);
%                         [updatedModel_PosRxn, updatedEColiKeggID, mainCouplingRxn_PosRxn] = restrictFluxOfAddedReaction(compound, updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, 1, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1});
                        [updatedModel_PosRxn, updatedEColiKeggID, mainCouplingRxn_PosRxn, posBoundVal] = setBoundsForAddedRxn(updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, 1, prodEnzymesList{prodIdx-2,1}, selectedOperators(operatorsPathwayDetailsIdx));
                        mutantPosFluxValue = optimizeCbModel(updatedModel_PosRxn, 'max', false, false);

                        % couple the reaction so the flux would be in the
                        % opposite direction as the original reaction flux
%                         [updatedModel_NegRxn, updatedEColiKeggID, mainCouplingRxn_NegRxn] = coupleRxnFlux(updatedModel, originalEColiKeggID, prodEnzymesList{prodIdx-2,1},...
%                                             pathwayFluxDetails, updatedRxnIDList(end), couplePercent, -1);

%                         [updatedModel_NegRxn, updatedEColiKeggID, mainCouplingRxn_NegRxn] = restrictFluxOfAddedReaction(compound, updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, -1, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1});
                        [updatedModel_NegRxn, updatedEColiKeggID, mainCouplingRxn_NegRxn, negBoundVal] = setBoundsForAddedRxn(updatedModel, originalEColiKeggID, updatedRxnIDList(end), couplePercent, -1, prodEnzymesList{prodIdx-2,1}, selectedOperators(operatorsPathwayDetailsIdx));
                        mutantNegFluxValue = optimizeCbModel(updatedModel_NegRxn, 'max', false, false);


                        % setting default for copuling to be in the same direction
                        coupleDirection = 1;
                        mutantFluxValue = mutantPosFluxValue; 
                        boundVal = posBoundVal;
                        
                        if (mutantPosFluxValue.f <= 0.0001 || mutantNegFluxValue.f <= 0.0001)
                            i;
                        end
                    
                        % check if the default coupling direction will change                
                        if mutantPosFluxValue.f < mutantNegFluxValue.f && mutantPosFluxValue.f > 0.0001
                            coupleDirection = 1;
                            mutantFluxValue = mutantPosFluxValue; 
                            boundVal = posBoundVal;
                        elseif mutantNegFluxValue.f < mutantPosFluxValue.f && mutantNegFluxValue.f > 0.0001
                            coupleDirection = -1;
                            mutantFluxValue = mutantNegFluxValue;
                            boundVal = negBoundVal;
                        end
                        modelAugmentFlag = false;
                    
                        if mutantFluxValue.f > 0.0001
                            index = find(EColiKeggID == prodCmpID(kfcProdIdx));
                            [rxnIDList, modelAugmentFlag] = updateMainModel(index, compound, prodCmpID(kfcProdIdx), rxnIDList, selectedOperators(operatorsPathwayDetailsIdx));
                        end

                        if modelAugmentFlag
                            stepsFBAResultsPos_OneStep = [stepsFBAResultsPos_OneStep; mutantPosFluxValue];
                            allModelPos_OneStep = [allModelPos_OneStep; updatedModel_PosRxn];

                            stepsFBAResultsNeg_OneStep = [stepsFBAResultsNeg_OneStep; mutantNegFluxValue];
                            allModelNeg_OneStep = [allModelNeg_OneStep; updatedModel_NegRxn];
                            
%                             [mainModel, EColiKeggID, mainCouplingRx] = restrictFluxOfAddedReaction(compound, mainModel, EColiKeggID, rxnIDList(end), couplePercent, coupleDirection, pathwayFluxDetails, prodEnzymesList{prodIdx-2,1});
                            [mainModel, EColiKeggID, mainCouplingRx] = setBoundsForAddedRxn_FinalModel(mainModel, EColiKeggID, rxnIDList(end), coupleDirection, prodEnzymesList{prodIdx-2,1}, boundVal, selectedOperators(operatorsPathwayDetailsIdx));
                            
%                             [mainModel, EColiKeggID, mainCouplingRxn] = coupleRxnFlux(mainModel, EColiKeggID,...
%                                 prodEnzymesList{prodIdx-2,1}, pathwayFluxDetails, rxnIDList(end), currentCouplePercent, coupleDirection); 
    %                         mutantFluxValue = optimizeCbModel(mainModel, 'max', false, false);
    %                         stepsFBAResults = [stepsFBAResults; mutantFluxValue];
    %                         originalRxn_couplingList = [originalRxn_couplingList; mainCouplingRxn];

                            selectedOperators(operatorsPathwayDetailsIdx).moleculePosition = products(prodIdx-2).ObtainedFromInput;
                            selectedOperators(operatorsPathwayDetailsIdx).coupleDirection = coupleDirection;
                            selectedOperators(operatorsPathwayDetailsIdx).OperatorIdx = operatorsPathwayDetailsIdx;

                            operatorsPerStepList = [operatorsPerStepList; selectedOperators(operatorsPathwayDetailsIdx)];
                            
                            allModel_Combined = [allModel_Combined; mainModel];

                            % saving the mol file in a different directory for
                            % later consideration
                            if kfcProdIdx == 1
                                newNameProdFile = [currentProdFile(1:end-4), '-', num2str(prodCmpID(kfcProdIdx)), '.mol'];
                                sourceFolder = [productsFolder, currentProdFile];
                                destinationFolder = [savedProductsFolder, newNameProdFile];
                                movefile(sourceFolder, destinationFolder);
                            end

                        end

                    end
                    
                     
                end
            end
%             end
            
        end
    end
end
end

function [phaseI, phaseImainAtom, phaseImainAtomPos] = sortRs(cypData)
% This function sorts the sturuct based on phase.Reactant.R
for i = 1:length(cypData)
    main{i} = cypData(i).Reactant.R;
end
[sortedR, order] = sort(main);
phaseI = cypData(order);
phaseImainAtom = unique(sortedR);
for i = 1:length(phaseImainAtom)
    for j = 1:length(sortedR)
        if (findstr(phaseImainAtom{i}, sortedR{j}))
            phaseImainAtomPos(i) = j;
            break;
        end
    end
end
end

% This function is to update the mainModel used with the product and
% reactions related to it. There are three different cases related to
% identified products, Case A, B and C explained below in appropriate spots

function [rxnIDList,modelAugmentFlag] = updateMainModel(productIdx, queryCompound, prodCmpID, rxnIDList, operatorsPathwayDetails)
global EColiKeggID
global mainModel
global substrateIDs
global prodIDsList

modelAugmentFlag = false;

cmpIdx = find(EColiKeggID == queryCompound);
if length(cmpIdx)>1
    cmpNames = mainModel.mets(cmpIdx);
    for i = 1:length(cmpNames)
        cmpStr = char(cmpNames(i,end)); 
        if  cmpStr(end) == 'c' 
            cmpIdxTemp = cmpIdx(i);
        end
    end
    cmpIdx = cmpIdxTemp;
end


if ~isempty(productIdx)
   % if the product already exists in ecoli check if there is a
   % reaction connecting the product and the current Ecoli compound
   % being queried
    if length(productIdx)>1
        prodNames = mainModel.mets(productIdx);
        for i = 1:length(prodNames)
            prodStr = char(prodNames(i,end)); 
            if  prodStr(end) == 'c' 
                prodIdxTemp = productIdx(i);
            end
        end
        productIdx = prodIdxTemp;
    end
    
    for i = 1:length(cmpIdx)
        cmpRxn = find(mainModel.S(cmpIdx(i),:)~=0);
        for j = 1:length(productIdx)
             prodRxn = find(mainModel.S(productIdx(j),:)~=0);
             
             % find reactions that are common between compound and product
             commonRxns = intersect(cmpRxn, prodRxn);
             
             if ~isempty(commonRxns)
                 [flag, rxn] = areCompoundsConnected(commonRxns, cmpIdx(i), productIdx(j));
                 
                 % Case A: product is in ecoli, and there is a rxn
                 % connecting the product and compound, so update the rxn
                 % capacity in the model
                 if flag
                     flag = checkSameCofactors(mainModel, rxn, operatorsPathwayDetails);
                     
                     if flag
                         mainModel.lb(rxn) = 1.1 * mainModel.lb(rxn);
                         mainModel.ub(rxn) = 1.1 * mainModel.ub(rxn);
                         rxnIDList = [rxnIDList; rxn]; 
                         substrateIDs = [substrateIDs; queryCompound];
                         prodIDsList = [prodIDsList; prodCmpID];
                         modelAugmentFlag = true;
                     else
                         rxnIDList = addRxnToModel(cmpIdx(i), productIdx(j), prodCmpID, rxnIDList, operatorsPathwayDetails);
                         substrateIDs = [substrateIDs; queryCompound];
                         prodIDsList = [prodIDsList; prodCmpID];
                         modelAugmentFlag = true;
                     end 
                 end
                 
             % Case B: product is in ecoli, but there is no rxn connecting
             % the product and compound, so update the model by adding
             % a rxn to connect product and compound
             else
                 rxnIDList = addRxnToModel(cmpIdx(i), productIdx(j), prodCmpID, rxnIDList, operatorsPathwayDetails);
                 substrateIDs = [substrateIDs; queryCompound];
                 prodIDsList = [prodIDsList; prodCmpID];
                 modelAugmentFlag = true;
             end
        end
    end

% Case C: product is not in ecoli, so add product to mets
% list and add a rxn connecting compound and product    
else
    if length(cmpIdx) > 1
        [prodIdx, rxnIDList] = addProduct_RxnToModel(cmpIdx(1), prodCmpID, rxnIDList, operatorsPathwayDetails);
        substrateIDs = [substrateIDs; queryCompound];
        prodIDsList = [prodIDsList; prodCmpID];
        EColiKeggID(end+1) = prodCmpID;
        for i = 2:length(cmpIdx)
            rxnIDList = addRxnToModel(cmpIdx(i), prodIdx, prodCmpID, rxnIDList, operatorsPathwayDetails);
            substrateIDs = [substrateIDs; queryCompound];
            prodIDsList = [prodIDsList; prodCmpID];
            modelAugmentFlag = true;
        end
        
    else
        EColiKeggID(end+1) = prodCmpID;
        [prodIdx, rxnIDList] = addProduct_RxnToModel(cmpIdx(1), prodCmpID, rxnIDList, operatorsPathwayDetails);
        substrateIDs = [substrateIDs; queryCompound];
        prodIDsList = [prodIDsList; prodCmpID];
        modelAugmentFlag = true;
    end
end
end

% This function is to update the mainModel used with the product and
% reactions related to it. There are three different cases related to
% identified products, Case A, B and C explained below in appropriate spots

function [updatedRxnIDList, mainModel, modelAugmentFlag] = updateModel_OneStep(origianlModel, originalEColiKeggID, productIdx, queryCompound, prodCmpID, selectedOperators)

mainModel = origianlModel;
updatedRxnIDList = [];
modelAugmentFlag = false;
cmpIdxTemp = [];

% make sure if the compound exists in different compartment to choose the
% one located in cytoplasem --> ending with "_c"
cmpIdx = find(originalEColiKeggID == queryCompound);
if length(cmpIdx)>1
    cmpNames = mainModel.mets(cmpIdx);
    for i = 1:length(cmpNames)
        cmpStr = char(cmpNames(i,end)); 
        if  cmpStr(end) == 'c' 
            cmpIdxTemp = cmpIdx(i);
        end
    end
    if (~isempty(cmpIdxTemp))
        cmpIdx = cmpIdxTemp;
    else
        return;
    end
end

modelAugmentFlag = true;
if ~isempty(productIdx)
    
   % if the product already exists in ecoli check if there is a
   % reaction connecting the product and the current Ecoli compound
   % being queried
    if length(productIdx)>1
        prodNames = mainModel.mets(productIdx);
        for i = 1:length(prodNames)
            prodStr = char(prodNames(i,end)); 
            if  prodStr(end) == 'c' 
                prodIdxTemp = productIdx(i);
            end
        end
        productIdx = prodIdxTemp;
    end
    
    for i = 1:length(cmpIdx)
        cmpRxn = find(mainModel.S(cmpIdx(i),:)~=0);
        for j = 1:length(productIdx)
             prodRxn = find(mainModel.S(productIdx(j),:)~=0);
             
             % find reactions that are common between compound and product
             commonRxns = intersect(cmpRxn, prodRxn);
             
             if ~isempty(commonRxns)
                 [flag, rxn] = areCompoundsConnected(commonRxns, cmpIdx(i), productIdx(j));
                 
                 % Case A: product is in ecoli, and there is a rxn
                 % connecting the product and compound, so update the rxn
                 % capacity in the model
                 if flag
                     
                     % check if the equation has the same cofactors as the
                     % original coupled equation
                     flag = checkSameCofactors(mainModel, rxn, selectedOperators);
                     
                     if flag
                         mainModel.lb(rxn) = 1.1 * mainModel.lb(rxn);
                         mainModel.ub(rxn) = 1.1 * mainModel.ub(rxn);
                         updatedRxnIDList = rxn;
                         
                     % Case B: product is in ecoli, but there is no rxn connecting
                     % the product and compound, so update the model by adding
                     % a rxn to connect product and compound
                     else
                         [nmets, nrxns] = size(mainModel.S);
                         mainModel.S(cmpIdx(i), nrxns+1) = -1;
                         mainModel.S(productIdx(j), nrxns+1) = 1;
                         
                         if ~isempty(selectedOperators)
                             mainModel = addCofactorsToModel(mainModel, originalEColiKeggID, selectedOperators, nrxns);
                         end
                         
                         productRxnStrName = int2str(prodCmpID);
                         for charIdx = 5:-1:length(productRxnStrName)
                             productRxnStrName = ['0', productRxnStrName];
                         end
                         mainModel.rxns{nrxns+1,1} = ['ExRxn_C', productRxnStrName];
                         mainModel.rxnNames{nrxns+1,1} = ['ExRxn_C', productRxnStrName];
                         mainModel.lb(nrxns+1) = -1000;
                         mainModel.ub(nrxns+1) = 1000;
                         mainModel.c(nrxns+1) = 0;
                         mainModel.rev(nrxns+1) = 1;
                         updatedRxnIDList = nrxns+1;
                         
                     end
                 end
                     
             % Case B: product is in ecoli, but there is no rxn connecting
             % the product and compound, so update the model by adding
             % a rxn to connect product and compound
             else
                 [nmets, nrxns] = size(mainModel.S);
                 mainModel.S(cmpIdx(i), nrxns+1) = -1;
                 mainModel.S(productIdx(j), nrxns+1) = 1;
                 
                if ~isempty(selectedOperators)
                    mainModel = addCofactorsToModel(mainModel, originalEColiKeggID, selectedOperators, nrxns);
                end
                 
                 productRxnStrName = int2str(prodCmpID);
                 for charIdx = 5:-1:length(productRxnStrName)
                     productRxnStrName = ['0', productRxnStrName];
                 end
                 mainModel.rxns{nrxns+1,1} = ['ExRxn_C', productRxnStrName];
                 mainModel.rxnNames{nrxns+1,1} = ['ExRxn_C', productRxnStrName];
                 mainModel.lb(nrxns+1) = -1000;
                 mainModel.ub(nrxns+1) = 1000;
                 mainModel.c(nrxns+1) = 0;
                 mainModel.rev(nrxns+1) = 1;
                 updatedRxnIDList = nrxns+1;
                
             end
        end
    end

% Case C: product is not in ecoli, so add product to mets
% list and add a rxn connecting compound and product    
else
    if length(cmpIdx) > 1
        
        cmpNames = model.mets(cmpIdx);
        
        % find the index of the compound that is in the cytoplacem compartment
        for j = 1:length(cmpNames)
            cmpStr = char(cmpNames(j,end));
            if  cmpStr(end) == 'c'
                cmpIdxTemp = cmpIdx(j);
                break;
            end
        end
        cmpIdx = cmpIdxTemp;
        
        [nmets, nrxns] = size(mainModel.S);
        mainModel.S(cmpIdx(1), nrxns+1) = -1;
        mainModel.S(nmets+1, nrxns+1) = 1;
        
        if ~isempty(selectedOperators)
            mainModel = addCofactorsToModel(mainModel, originalEColiKeggID, selectedOperators, nrxns);
        end
        
        productRxnStrName = int2str(prodCmpID);
        for charIdx = 5:-1:length(productRxnStrName)
            productRxnStrName = ['0', productRxnStrName];
        end
        mainModel.rxns{nrxns+1,1} = ['Rxn_C', productRxnStrName];
        mainModel.rxnNames{nrxns+1,1} = ['Rxn_C', productRxnStrName];
        mainModel.mets{nmets+1} = ['C', productRxnStrName];
        mainModel.lb(nrxns+1) = -1000;
        mainModel.ub(nrxns+1) = 1000;
        mainModel.c(nrxns+1) = 0;
        mainModel.rev(nrxns+1) = 1;
        mainModel.b(nmets+1) = 0;
        mainModel.csense = [mainModel.csense 'E'];
        updatedRxnIDList = nrxns+1;
        
        % adding an outgoing/Exchange reaction for the new products that are added
        % to the model to make sure if they are produced then the flux will have a
        % way out of the cell
        mainModel.S(nmets+1, nrxns+2) = -1;
        mainModel.rxns{nrxns+2,1} = ['ExRxn_C', productRxnStrName];
        mainModel.rxnNames{nrxns+2,1} = ['ExRxn_C', productRxnStrName];
        mainModel.lb(nrxns+2) = 0;
        mainModel.ub(nrxns+2) = 1000;
        mainModel.c(nrxns+2) = 0;
        mainModel.rev(nrxns+2) = 0;

%         for i = 2:length(cmpIdx)
%             [nmets, nrxns] = size(mainModel.S);
%              mainModel.S(cmpIdx(i), nrxns+1) = -1;
%              mainModel.S(prodIdx, nrxns+1) = 1;
%                 
%              productRxnStrName = int2str(prodCmpID);
%              for charIdx = 5:-1:length(productRxnStrName)
%                  productRxnStrName = ['0', productRxnStrName];
%              end
%              mainModel.rxns{nrxns+1,1} = ['ExRxn_C', productRxnStrName];
%              mainModel.rxnNames{nrxns+1,1} = ['ExRxn_C', productRxnStrName];
%              mainModel.lb(nrxns+1) = -1000;
%              mainModel.ub(nrxns+1) = 1000;
%              mainModel.c(nrxns+1) = 0;
%              mainModel.rev(nrxns+1) = 1;
% %              updatedRxnIDList = nrxns+1;
%              
%         end
        
    else
        
        [nmets, nrxns] = size(mainModel.S);
        mainModel.S(cmpIdx(1), nrxns+1) = -1;
        mainModel.S(nmets+1, nrxns+1) = 1;
        
        if ~isempty(selectedOperators)
            mainModel = addCofactorsToModel(mainModel, originalEColiKeggID, selectedOperators, nrxns);
        end
                 
        productRxnStrName = int2str(prodCmpID);
        for charIdx = 5:-1:length(productRxnStrName)
            productRxnStrName = ['0', productRxnStrName];
        end
        mainModel.rxns{nrxns+1,1} = ['Rxn_C', productRxnStrName];
        mainModel.rxnNames{nrxns+1,1} = ['Rxn_C', productRxnStrName];
        mainModel.mets{nmets+1} = ['C', productRxnStrName];
        mainModel.lb(nrxns+1) = -1000;
        mainModel.ub(nrxns+1) = 1000;
        mainModel.c(nrxns+1) = 0;
        mainModel.rev(nrxns+1) = 1;
        mainModel.b(nmets+1) = 0;
        mainModel.csense = [mainModel.csense 'E'];
        updatedRxnIDList = nrxns+1;
        
        % adding an outgoing/Exchange reaction for the new products that are added
        % to the model to make sure if they are produced then the flux will have a
        % way out of the cell
        mainModel.S(nmets+1, nrxns+2) = -1;
        mainModel.rxns{nrxns+2,1} = ['ExRxn_C', productRxnStrName];
        mainModel.rxnNames{nrxns+2,1} = ['ExRxn_C', productRxnStrName];
        mainModel.lb(nrxns+2) = 0;
        mainModel.ub(nrxns+2) = 1000;
        mainModel.c(nrxns+2) = 0;
        mainModel.rev(nrxns+2) = 0;

    end
end
end

function [flag, rxn] = areCompoundsConnected(commonRxns, cmpIdx, productIdx)
global mainModel
flag = false;
rxn = [];


for i = 1:length(commonRxns)
    [~, ~, compCoeff] = find(mainModel.S(cmpIdx, commonRxns(i)));
    [~, ~, prodCoeff] = find(mainModel.S(productIdx, commonRxns(i)));
    if compCoeff == 1 && prodCoeff == -1
        flag = true;
        rxn = commonRxns(i);
        return;
        
    elseif compCoeff == -1 && prodCoeff == 1
        flag = true;
        rxn = commonRxns(i);
        return;
        
    end
end
end

function rxnIDList = addRxnToModel(cmpIdx, productIdx, prodCmpID, rxnIDList, operatorsPathwayDetails)
global mainModel
global EColiKeggID

[nmets, nrxns] = size(mainModel.S);
mainModel.S(cmpIdx, nrxns+1) = -1;
mainModel.S(productIdx, nrxns+1) = 1;

if ~isempty(operatorsPathwayDetails)
    mainModel = addCofactorsToModel(mainModel, EColiKeggID, operatorsPathwayDetails, nrxns);
end

productRxnStrName = int2str(prodCmpID);
for i = 5:-1:length(productRxnStrName)
    productRxnStrName = ['0', productRxnStrName];
end
mainModel.rxns{nrxns+1,1} = ['Rxn_C', productRxnStrName, '_R', int2str(nrxns+1)];
mainModel.rxnNames{nrxns+1,1} = ['Rxn_C', productRxnStrName, '_R', int2str(nrxns+1)];
mainModel.lb(nrxns+1) = -1000;
mainModel.ub(nrxns+1) = 1000;
mainModel.c(nrxns+1) = 0;
mainModel.rev(nrxns+1) = 1;
rxnIDList = [rxnIDList; nrxns+1];

end

function [nmets, rxnIDList] = addProduct_RxnToModel(cmpIdx, prodCmpID, rxnIDList, operatorsPathwayDetails)
global mainModel
global EColiKeggID
 
[nmets, nrxns] = size(mainModel.S);
rxnIDList = [rxnIDList; nrxns+1];
mainModel.S(cmpIdx, nrxns+1) = -1;
mainModel.S(nmets+1, nrxns+1) = 1;


if ~isempty(operatorsPathwayDetails)
    mainModel = addCofactorsToModel(mainModel, EColiKeggID, operatorsPathwayDetails, nrxns);
end

productRxnStrName = int2str(prodCmpID);
for i = 5:-1:length(productRxnStrName)
    productRxnStrName = ['0', productRxnStrName];
end
mainModel.rxns{nrxns+1,1} = ['Rxn_C', productRxnStrName, '_R', int2str(nrxns+1)];
mainModel.rxnNames{nrxns+1,1} = ['Rxn_C', productRxnStrName, '_R', int2str(nrxns+1)];
mainModel.mets{nmets+1} = ['C', productRxnStrName];
mainModel.lb(nrxns+1) = -1000;
mainModel.ub(nrxns+1) = 1000;
mainModel.c(nrxns+1) = 0;
mainModel.rev(nrxns+1) = 1;
mainModel.b(nmets+1) = 0;
mainModel.csense = [mainModel.csense 'E'];

% adding an outgoing/Exchange reaction for the new products that are added
% to the model to make sure if they are produced then the flux will have a
% way out of the cell
mainModel.S(nmets+1, nrxns+2) = -1;
mainModel.rxns{nrxns+2,1} = ['ExRxn_C', productRxnStrName];
mainModel.rxnNames{nrxns+2,1} = ['ExRxn_C', productRxnStrName];
mainModel.lb(nrxns+2) = 0;
mainModel.ub(nrxns+2) = 1000;
mainModel.c(nrxns+2) = 0;
mainModel.rev(nrxns+2) = 0;

nmets = nmets+1;
end