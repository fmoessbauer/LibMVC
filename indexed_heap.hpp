#ifndef INDEXED_HEAP_INCLUDED
#define INDEXED_HEAP_INCLUDED

#include <vector>
#include <map>
#include <functional>

template<
  typename T,
  typename Container = std::vector<T>,
  typename Index     = std::map<T, typename Container::size_type>,
  class    Compare   = std::less<T>>
class Indexed_Heap {
  public:
    using container_type  = Container;
    using value_compare   = Compare;
    using value_type      = typename Container::value_type;
    using size_type       = typename Container::size_type;
    using reference       = typename Container::reference;
    using const_reference = typename Container::const_reference;

  private:
    Container data;
    Index     index;
    Compare   comp;

  private:
    void element_swap(
        const size_type & a,
        const size_type & b)
    {
      std::swap(data[a], data[b]);
      index.at(data[a]) = a;
      index.at(data[b]) = b;
    }

    bool is_leaf(const size_type & pos) const {
      if( (pos >= data.size()/2) && (pos<data.size()) ) return true;
      return false;
    }

    size_type left_child(const size_type & pos) const {
      return 2*pos+1;
    }

    size_type right_child(const size_type & pos) const {
      return 2*pos+2;
    }

    size_type parent(const size_type & pos) const {
      return (pos-1)/2;
    }

    void shift_down(size_type pos){
      while(!is_leaf(pos)){
        auto j = left_child(pos);
        auto rc = right_child(pos);
        if( (rc<data.size()) && (comp(data[rc], data[j])) ) j = rc;
        if( comp(pos, j)) return;
        element_swap(pos, j);
        pos = j;
      }
    }

  public:
    explicit Indexed_Heap(const Compare & compare = Compare())
    : comp(compare){ }

  public:
    const_reference top() const {
      return data[0];
    }

    bool empty() const {
      return data.size() == 0;
    }

    size_type size() const {
      return data.size();
    }

    void push(const value_type & value){
      int curr     = data.size();
      index[value] = curr;
      data.push_back(value);

      while(curr!=0 && comp(data[curr], data[parent(curr)])){
        element_swap(curr, parent(curr));
        curr = parent(curr);
      }
    }

    void pop(){
      element_swap(0, data.size()-1);
      index.erase(data.back());
      data.pop_back();
      if(data.size() != 0) shift_down(0);
    }

    /**
     * returns 1 if element is in heap, 0 otherwise
     */
    size_type count(const value_type & value) const {
      return index.count(value);
    }

    /**
     * remove element at given position and restore heap
     */
    void erase(const size_type & pos){
      if(pos == (data.size()-1)){
        index.erase(data.back());
        data.pop_back();
      } else {
        element_swap(pos, data.size()-1);
        index.erase(data.back());
        data.pop_back();
        auto idx = pos;
        while( (idx !=0 ) && (comp(data[idx],data[parent(idx)]))){
          element_swap(idx, parent(idx));
          idx = parent(idx);
        }
        if(data.size()!=0) shift_down(idx);
      }
    }

    /**
     * return position in heap of given value
     */
    size_type operator[](const value_type & value) const {
      return index.at(value);
    }
};
#endif

