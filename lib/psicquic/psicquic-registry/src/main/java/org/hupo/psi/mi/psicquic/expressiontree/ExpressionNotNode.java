package org.hupo.psi.mi.psicquic.expressiontree;

import org.hupo.psi.mi.psicquic.ols.client.SelfDiscoveringOntologyTree;
import org.hupo.psi.mi.psicquic.registry.ServiceType;
/**
 * 
 */

/**
 * @author Erik Pfeiffenberger
 *
 */
public class ExpressionNotNode implements ExpressionNode{
	
	private ExpressionNode child;
	
	public ExpressionNotNode(ExpressionNode child) {
		this.child = child;
	}

	public boolean evaluate(ServiceType service, SelfDiscoveringOntologyTree sdoTree) {
		return !child.evaluate(service, sdoTree);
	}

	public ExpressionNode getChild() {
		return child;
	}

	public void setChild(ExpressionNode child) {
		this.child = child;
	}

	public String getRepresentation() {
		return "NOT";
	}

	public ExpressionNode getLeftChild() {
		return getChild();
	}

	public ExpressionNode getRightChild() {
		return null;
	}


}
