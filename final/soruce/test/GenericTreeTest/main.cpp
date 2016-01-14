#include <bits/stdc++.h>
using namespace std;

struct TreeData {
	uint32_t v[4];
};
struct TreeNode {
	// node link information
	TreeNode *parent, *childFirst, *childTail;
	TreeNode *siblingNext, *siblingPrev;
	// node data
	TreeData e;
};
struct LinearTreeNode {
	// node link information
	uint32_t parentOffset;
	uint32_t childOffsetHead, childOffsetTail;
	uint32_t siblingOffsetNext, siblingOffsetPrev;
	// faster record
	uint32_t numChild;
	// node data
	TreeData e;
	uint32_t visitCount;
};

int nodesize = 0;
TreeNode* parser() {
	nodesize++;
	TreeNode *u = (TreeNode *) malloc(sizeof(TreeNode));
	u->parent = u->childFirst = u->childTail = NULL;
	u->siblingNext = u->siblingPrev = NULL;
	
	TreeNode *p = NULL;
	int has = 0;
	for (; ; ) {
		scanf("%d", &has);
		if (!has)
			break;
		TreeNode *q = parser();
		if (p == NULL) {
			u->childFirst = q;
		} else {
			p->siblingNext = q, q->siblingPrev = p;
		}
		
		p = q;
	}
	u->childTail = p;
	return u;
}
uint32_t flattenTree(TreeNode *node, uint32_t *offset, uint32_t parentOffset, LinearTreeNode *_mem) {
	LinearTreeNode *linearNode = &_mem[*offset];
	uint32_t myOffset = (*offset)++;
	linearNode->childOffsetHead = linearNode->childOffsetTail = -1;
	linearNode->parentOffset = parentOffset;
	linearNode->siblingOffsetNext = linearNode->siblingOffsetPrev = -1;
	if (node->childFirst != NULL) {
		TreeNode *ch = node->childFirst;
		uint32_t prevch = -1;
		for (; ch != NULL; ch = ch->siblingNext) {
			uint32_t nextch = flattenTree(ch, offset, myOffset, _mem);
			if (prevch != -1) {	
				LinearTreeNode *p = &_mem[prevch];
				LinearTreeNode *q = &_mem[nextch];
				p->siblingOffsetNext = nextch;
				q->siblingOffsetPrev = prevch;
			} else {
				linearNode->childOffsetHead = nextch;
			}
			prevch = nextch;
		}
		linearNode->childOffsetTail = prevch;
	}
	return myOffset;
}
uint32_t hash1 = 0, hash2 = 0, order = 0;
void recursiveTraversal(LinearTreeNode *node, LinearTreeNode *_mem) {
	hash1 += *(uint32_t *) node * (order++);
#ifdef DEBUG
	printf("%p\n", node);
#endif
	uint32_t offset = node->childOffsetHead;
	if (offset == -1)
		return ;
	for (LinearTreeNode *u; offset != -1; offset = u->siblingOffsetNext) {
		u = &_mem[offset];
#ifdef DEBUG
//		printf("%p %p\n", node, u);
#endif
		recursiveTraversal(u, _mem);
	}
}
void iteratorTraversal(uint32_t offset, LinearTreeNode *_mem) {
	bool process = true;
	while (offset != -1) {
		LinearTreeNode *node = &_mem[offset];
		if (process) {
#ifdef DEBUG
			printf("%p\n", node);
#endif
			hash2 += *(uint32_t *) node * (order++);
		}
		if (node->childOffsetHead != -1 && process) {
			offset = node->childOffsetHead;
			process = true;
		} else if (node->siblingOffsetNext != -1) {
			offset = node->siblingOffsetNext;
			process = true;
		} else {
			offset = node->parentOffset;
			process = false;
		}
	}
}
int main() {
	nodesize = 0;
	scanf("%*d");
	TreeNode *root = parser();
	LinearTreeNode *_mem = (LinearTreeNode *) malloc(sizeof(LinearTreeNode) * nodesize);
	
	uint32_t offset = 0;
	flattenTree(root, &offset, -1, _mem);
	hash1 = hash2 = 0;
	printf("sizeof(LinearTreeNode) = %d\n", sizeof(LinearTreeNode));
#define MAXLOOP 10000
#ifdef RECTEST
	for (int i = 0; i < MAXLOOP; i++)
		recursiveTraversal(&_mem[0], _mem);
	printf("%lu\n", hash1);
#endif
#ifdef ITETEST
	for (int i = 0; i < MAXLOOP; i++)
		iteratorTraversal(0, _mem);
	printf("%lu\n", hash2);
#endif
	return 0;
}
/*
1
	1 
		1 0 
		1 0 
		1 0 
		0
	1 
		1 0 
		1 0 
		0 
	0
*/
