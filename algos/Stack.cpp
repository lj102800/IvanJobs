/**
 * demo for stack
 */
#define MAX 100
int c[MAX];

// index starting from 1, and we will not use index 0.
class Stack {
public:
	Stack(int size);
	bool push(int v);
	bool pop();
	bool isEmpty();
	bool isFull();
private:
	int top;
	int size;
};

Stack::Stack(int size) {
	this->size = size;
	this->top = 0;
}

bool Stack::push(int v) {
	if (this->isFull())
		return false;
	this->top++;
	c[this->top] = v;
	return true;
}

bool Stack::pop() {
	if (this->isEmpty())
		return false;
	this->top--;
}

bool Stack::isEmpty() {
	if (this->top == 0)
		return true;
	else 
		return false;
}

bool Stack::isFull() {
	if (this->top == this->size)
		return true;
	else 
		return false;
}