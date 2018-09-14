% Function that returns all enzymes catalyzing each reaction in a rxnsList
% input:
%   rxns: a list of the model's KEGG reactions IDs
% output: 
%   enzymeECNumList: a list of EC numbers catalyzing each reaction
function enzymeECNumList = getEnzymeIDInEcoli(rxns)
global reactions
global indexReactions

enzymeECNumList = [];

counter = 1;
for j = 1: length(rxns)
    % for each reaction in the pathway, find all EC numbers of enzymes
    % catalyzing the reaction
    if rxns(j) > 0
        reactionID = rxns(j);
        reactionIdx = indexReactions(reactionID);
        enzymeECNumList{j}{counter} = reactions(reactionIdx).Enzymes;
    end
end

end