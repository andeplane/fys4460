#pragma once
#include <cstddef>
#include <iostream>   // one of two
 using namespace std;    // two of two, yay!

template <class T> class Node
{
public:
    T data; //the object information
    Node<T>* next; //pointer to the next node element
    Node<T>* prev; //pointer to the previous node element

public:
    Node() {
        next = NULL;
        prev = NULL;
    }

    //Methods omitted for brevity
};

template <class T> class LinkedList
{
    Node<T> *head, *tail;
public:
    int count;

    LinkedList()
    {
        head = NULL;
        tail = NULL;
        count = 0;
    }

    void add(T elem) {
        if(head == NULL) {
            head = new Node<T>;
            head->data = elem;
            tail = head;
        }
        else {
            Node<T> *temp = new Node<T>;
            temp->data = elem;

            tail->next = temp; // Set this element as tail's next-element
            temp->prev = tail; // Set tail as element's prev-element

            tail = temp; // Update tail
        }
        count++;
    }

    void remove(Node<T>* node) {
        if(node == tail) {
            tail = tail->prev;
            tail->next = NULL;
        } else if(node == head) {
            head = head->next;
            head->prev = NULL;
        }
        else {
            node->next->prev = node->prev;
            node->prev->next = node->next;
        }
        delete node;

        count--;
    }

    void remove(T elem) {
        Node<T> *node = head;

        while(node != NULL) {
            if(node->data == elem) {
                if(node == head) head = node->next;
                if(node == tail) tail = node->prev;

                if(node->prev != NULL) node->prev->next = node->next;
                if(node->next != NULL) node->next->prev = node->prev;
                count--;
                delete node;

                return;
            }

            node = node->next;
        }
    }

};
