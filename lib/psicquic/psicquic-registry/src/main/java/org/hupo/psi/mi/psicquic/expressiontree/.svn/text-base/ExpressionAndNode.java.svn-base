package org.hupo.psi.mi.psicquic.expressiontree;

import org.hupo.psi.mi.psicquic.ols.client.SelfDiscoveringOntologyTree;
import org.hupo.psi.mi.psicquic.registry.ServiceType;

public class ExpressionAndNode implements ExpressionNode{
	private ExpressionNode leftChild;
	private ExpressionNode rightChild;
	
	public ExpressionAndNode(ExpressionNode leftChild, ExpressionNode rightChild) {
		this.leftChild = leftChild;
		this.rightChild = rightChild;
	}

	public boolean evaluate(ServiceType service, SelfDiscoveringOntologyTree sdoTree) {
		return leftChild.evaluate(service, sdoTree) && rightChild.evaluate(service, sdoTree);
	}

	public String getRepresentation() {
		return "AND";
	}

	public ExpressionNode getLeftChild() {
		return leftChild;
	}

	public ExpressionNode getRightChild() {
		return rightChild;
	}


}
