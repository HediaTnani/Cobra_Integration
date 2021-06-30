
initCobraToolbox()
cd ~/cobra_integration/
fileName = 'iMM1865.xml';
model = readCbModel(fileName);
expr=readtable('MeanExpressionValues_entrez.csv')
expressionData.gene=table2array(expr(:,2))
exp=table2array(expr(:,3))
expressionData.value=exp;
[expressionRxns, parsedGPR] = mapExpressionToReactions(model, expressionData)
