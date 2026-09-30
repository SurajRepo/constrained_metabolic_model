function mapExpressionsToReactions(model, genes, values)
    print("Mapping transcriptomics data to reactions \n" )
    gene_reaction = []
    for reaction in reactions(model)
        orcurExpr= [];
        associated_genes = model.mat["grRules"][findfirst(==(reaction), reactions(model))]
        if !isempty(associated_genes)
            curExprArr = COBREXA._parse_grr(associated_genes)
            for expR in curExprArr
                if length(expR)>=1
                    #And condition
                    andCurExp = []
                    for innerExpr in expR
                        geneIndex= findfirst(==(innerExpr), genes)
                        if !isnothing(geneIndex)
                            push!(andCurExp,values[geneIndex])
                        end
                    end
                    if !isempty(andCurExp)
                        push!(orcurExpr,minimum(andCurExp))  
                    end            
                end
            end
        end
        #or condition
        if isempty(orcurExpr)
            push!(gene_reaction,1)
        else
            push!(gene_reaction,maximum(orcurExpr))
        end  
    end
    return gene_reaction
end

function constrainReactions(model, geneReactionExpression, gamma)
    reaction =  reactions(model)
    bound = bounds(model)
    default_lower_bound = bound[1]
    default_upper_bound = bound[2]
    model = convert(CoreModel, model)
    for index in 1:length(geneReactionExpression)
        if geneReactionExpression[index] !=0
            if geneReactionExpression[index] >= 1
                lower_bound = default_lower_bound[index]*(1 + gamma*log(geneReactionExpression[index]))
                upper_bound = default_upper_bound[index]*(1 + gamma*log(geneReactionExpression[index]))
            else
                lower_bound = default_lower_bound[index]/(1+gamma*abs(log(geneReactionExpression[index])));
                upper_bound = default_upper_bound[index]/(1+gamma*abs(log(geneReactionExpression[index])));        
            end    
            #print(lower_bound, '\n')
            #print(reaction[index], '\t', default_upper_bound[index], '\t', geneReactionExpression[index], '\t', (1 + gamma*log(geneReactionExpression[index])), '\t', upper_bound , '\n')
            model = change_bound(model, string(reaction[index]), lower=lower_bound, upper=upper_bound)  
        end      
    end
    print("Finished constraining model with transcriptomics data \n" )
    return model
end

