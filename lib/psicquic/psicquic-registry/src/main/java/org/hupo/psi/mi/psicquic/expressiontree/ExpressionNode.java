package org.hupo.psi.mi.psicquic.expressiontree;

import org.hupo.psi.mi.psicquic.ols.client.SelfDiscoveringOntologyTree;
import org.hupo.psi.mi.psicquic.registry.ServiceType;

public interface ExpressionNode {

    public boolean evaluate(ServiceType service, SelfDiscoveringOntologyTree sdoTree);

    public String getRepresentation();
	
    public ExpressionNode getLeftChild();
    
	public ExpressionNode getRightChild();

}
