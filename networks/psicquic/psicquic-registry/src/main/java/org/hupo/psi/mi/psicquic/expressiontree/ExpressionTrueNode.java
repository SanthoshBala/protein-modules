/**
 * 
 */
package org.hupo.psi.mi.psicquic.expressiontree;

import org.hupo.psi.mi.psicquic.ols.client.SelfDiscoveringOntologyTree;
import org.hupo.psi.mi.psicquic.registry.ServiceType;

/**
 * @author Erik Pfeiffenberger
 *
 */
public class ExpressionTrueNode implements ExpressionNode{
	public ExpressionTrueNode() {
		// does nothing...
	}

	public boolean evaluate(ServiceType service, SelfDiscoveringOntologyTree sdoTree) {
		return true;
	}

	public String getRepresentation() {
		return "TRUE";
	}

	public ExpressionNode getLeftChild() {
		
		return null;
	}

	public ExpressionNode getRightChild() {
		return null;
	}

}
