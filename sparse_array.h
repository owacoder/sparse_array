/*
 * sparse_array.h
 *
 * Copyright © 2018 Oliver Adams
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the “Software”), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <typeinfo>

template<typename Key, typename T, typename DiffType, typename Generator> class sparse_array;

namespace impl {
    template<typename Key, typename T, typename DiffType, typename Generator>
    class GeneratorHelper
    {
        T val;
        Generator generator;

    public:
        GeneratorHelper(const sparse_array<Key, T, DiffType, Generator> *, const Generator *g) : val{}, generator(*g) {}

        const T &operator()(Key idx) {val = generator(idx); return val;}
    };

    template<typename Key, typename T, typename DiffType>
    class GeneratorHelper<Key, T, DiffType, void>
    {
        const sparse_array<Key, T, DiffType, void> *parent;

    public:
        GeneratorHelper(const sparse_array<Key, T, DiffType, void> *parent, const void *) : parent(parent) {}

        const T &operator()(Key);
    };

    class exception : public std::exception {
    public:
        exception(const char *what) : m(what) {}

        const char *what() const noexcept {return m;}

    private:
        const char *m;
    };
}

/*
 *  Key must have operator<() and operator==() defined, as well as the addition and subtraction operators
 *  Key also must be castable to DiffType and size_t (casting to size_t will only ever occur if <value of key> <= SIZE_MAX)
 *
 *  T may be of any type, as long as it is copy-constructable, default-constructable, and assignable
 *
 *  An invariant held throughout this implementation is that each map entry has at least one value stored in the
 *  "bucket", or vector, there.
 *
 *  For (very) sparse arrays (those without contiguous sections, or very few of them and small), the operations are:
 *
 *      Random element retrieval: O(logn)
 *      Random element existence: O(logn)
 *      Random element assigning: O(logn)
 *      Insertion at beginning: O(logn)
 *      Insertion in middle: O(logn)
 *      Insertion at end: O(logn)
 *      Erasing at beginning: O(logn)
 *      Erasing in middle: O(logn)
 *      Erasing at end: O(logn)
 *      Iteration: slower
 *      Everything else: O(1)
 *
 *  For very packed sparse arrays (those with massive contiguous sections), the operations are:
 *
 *      Random element retrieval: approaches O(1)
 *      Random element existence: approaches O(1)
 *      Random element assigning: approaches O(1)
 *      Insertion at beginning: approaches O(n)
 *      Insertion in middle: approaches O(n)
 *      Insertion at end: approaches O(1)
 *      Erasing at beginning: approaches O(n)
 *      Erasing in middle: approaches O(n)
 *      Erasing at end: approaches O(1)
 *      Iteration: faster
 *      Everything else: O(1)
 *
 * As you can see, read performance improves and insert and erasure performance declines
 * as contiguous sections grow, with the exception of inserting and deleting at the end.
 * Write performance for existing elements also improves when contiguous sections grow.
 *
 * All elements are strictly ordered, allowing O(n) ascending traversal of array keys.
 *
 */
template<typename Key, typename T, typename DiffType = uintmax_t, typename Generator = void>
class sparse_array {
    friend class const_iterator;

    static_assert(std::is_unsigned<DiffType>::value | std::is_class<DiffType>::value,
                  "sparse_array::DiffType must be an unsigned type or a user-defined class");

    size_t mElements;
    std::map<Key, std::vector<T>> mMap;
    mutable impl::GeneratorHelper<Key, T, DiffType, Generator> mGenerator;
    bool mUseRangeIters;
    T mDefaultValue;

    typedef std::map<Key, std::vector<T>> map_type;
    typedef typename map_type::iterator map_iterator;
    typedef typename map_type::const_iterator map_const_iterator;

public:
    class const_iterator : std::iterator<std::bidirectional_iterator_tag, T, std::ptrdiff_t, const T *, const T &>
    {
        friend class sparse_array;

        const sparse_array *mParent; // sparse array instance
        map_const_iterator mMapIter; // current map entry
        Key mKey; // Currently referenced key. May be invalid and out-of-bounds

        bool mIsRangeIter; // Set to true if non-existent (default-valued) elements should be iterated over

        const_iterator(const sparse_array *parent, map_const_iterator it, bool encompassAllValues)
            : mParent(parent)
            , mMapIter(it)
            , mKey{}
            , mIsRangeIter(encompassAllValues)
        {
            if (it != mParent->mMap.end())
                mKey = mMapIter->first;
        }

        const_iterator(const sparse_array *parent, map_const_iterator it, Key key, bool encompassAllValues)
            : mParent(parent)
            , mMapIter(it)
            , mKey(key)
            , mIsRangeIter(encompassAllValues)
        {}

        void advance() {
            if (mParent == nullptr)
                return;

            if (mIsRangeIter) {
                // If not at end, then check if we need to increment the map entry pointer
                // Otherwise, if at end, just increment
                if (mMapIter != mParent->mMap.end())
                {
                    auto nextIter = mMapIter;
                    ++nextIter;

                    // If the next entry is the end, bump mMapIter to the end
                    // Increment the key no matter what
                    if (nextIter == mParent->mMap.end()) {
                        if (mKey == mParent->vector_last_used_end_idx(mMapIter))
                            mMapIter = nextIter;

                        mKey = mKey + 1;
                    }
                    // Otherwise, are we ready to reset and move to next vector?
                    // If not, just increment
                    else {
                        mKey = mKey + 1;
                        if (mKey == nextIter->first)
                            mKey = (++mMapIter)->first;
                    }
                }
                else
                    mKey = mKey + 1;
            }
            else
            {
                // If at end of map entry, then obtain next map entry's initial value
                // If there is no next map entry, just increment
                if (mMapIter != mParent->mMap.end() && mKey == mParent->vector_last_used_end_idx(mMapIter) && ++mMapIter != mParent->mMap.end())
                    mKey = mMapIter->first;
                else
                    mKey = mKey + 1;
            }
        }

        void retreat() {
            if (mParent == nullptr)
                return;

            // If map entry is at end, then we need to reset back to the last map entry's final value
            // If there is no map entry, just decrement
            if (mMapIter == mParent->mMap.end()) {
                if (mMapIter != mParent->mMap.begin())
                    mKey = mParent->vector_last_used_end_idx(--mMapIter);
                else
                    mKey = mKey - 1;
            }
            // If key is equal to initial value of map entry, attempt to obtain previous map entry
            // Otherwise, if at beginning map entry, or not at beginning of current map entry, just decrement
            else if (mKey == mMapIter->first && mMapIter != mParent->mMap.begin()) {
                --mMapIter;

                // If a ranged iterator, just decrement,
                // Otherwise, get the last used index of the new map entry
                if (mIsRangeIter)
                    mKey = mKey - 1;
                else
                    mKey = mParent->vector_last_used_end_idx(mMapIter);
            }
            else
                mKey = mKey - 1;
        }

        bool atEnd() const {
            return mParent == nullptr || mMapIter == mParent->mMap.end();
        }

    public:
        typedef typename std::iterator<std::bidirectional_iterator_tag, T, std::ptrdiff_t, const T *, const T &>::reference reference;
        typedef typename std::iterator<std::bidirectional_iterator_tag, T, std::ptrdiff_t, const T *, const T &>::pointer pointer;

        const_iterator() : mParent(nullptr), mKey{} {}

        const_iterator &operator++() {advance(); return *this;}
        const_iterator &operator--() {retreat(); return *this;}
        const_iterator operator++(int) {const_iterator temp(*this); advance(); return temp;}
        const_iterator operator--(int) {const_iterator temp(*this); retreat(); return temp;}

        reference operator*() const
        {
            if (mIsRangeIter && mParent->is_not_in_vector(mMapIter, mKey))
                return mParent->mGenerator(mKey);

            return mMapIter->second[static_cast<size_t>(mKey - mMapIter->first)];
        }
        pointer operator->() const
        {
            if (mIsRangeIter && mParent->is_not_in_vector(mMapIter, mKey))
                return mParent->mGenerator(mKey);

            return std::addressof(mMapIter->second[static_cast<size_t>(mKey - mMapIter->first)]);
        }

        bool operator==(const const_iterator &other) const {
            if (atEnd() && other.atEnd())
                return true;

            return mParent == other.mParent &&
                    mMapIter == other.mMapIter &&
                    mKey == other.mKey;
        }

        bool operator!=(const const_iterator &other) const {
            return !(*this == other);
        }

        bool element_does_not_exist() const {return atEnd() || (mIsRangeIter && mParent->is_not_in_vector(mMapIter, mKey));}
        bool element_exists() const {return !element_does_not_exist();}
        const Key &index() const {return mKey;}
        reference value() const {return **this;}
    };

    typedef T value_type;

    sparse_array(const T &defaultValue, const std::vector<T> &args, bool iteratorsEncompassAllValuesInSpan = true)
        : mElements(0)
        , mGenerator(this, typeid(void) != typeid(Generator)? (throw impl::exception("sparse_array with non-default generator initialized without generator functor specified"), nullptr): nullptr)
        , mUseRangeIters(iteratorsEncompassAllValuesInSpan)
        , mDefaultValue(defaultValue)
    {
        mMap.insert(std::make_pair(0, args));
    }

    template<typename G>
    sparse_array(const std::vector<T> &args, const G &generator, bool iteratorsEncompassAllValuesInSpan = true)
        : mElements(0)
        , mGenerator(this, &static_cast<const Generator &>(generator))
        , mUseRangeIters(iteratorsEncompassAllValuesInSpan)
        , mDefaultValue{}
    {
        mMap.insert(std::make_pair(0, args));
    }

    sparse_array(const T &defaultValue, bool iteratorsEncompassAllValuesInSpan = true)
        : mElements(0)
        , mGenerator(this, typeid(void) != typeid(Generator)? (throw impl::exception("sparse_array with non-default generator initialized without generator functor specified"), nullptr): nullptr)
        , mUseRangeIters(iteratorsEncompassAllValuesInSpan)
        , mDefaultValue(defaultValue)
    {}

    template<typename G>
    sparse_array(const G &generator, bool iteratorsEncompassAllValuesInSpan = true)
        : mElements(0)
        , mGenerator(this, &static_cast<const Generator &>(generator))
        , mUseRangeIters(iteratorsEncompassAllValuesInSpan)
        , mDefaultValue{}
    {}

    // Complexity: O(1), constant-time
    // Set whether to use contiguous iterators (contiguous iterators include non-existent entries that take on the default value) for this class instance
    bool contiguous_iterators() const {return mUseRangeIters;}
    void set_contiguous_iterators(bool useContiguous) {mUseRangeIters = useContiguous;}

    // Returns either contiguous_xxx() or skip_xxx(), depending on the current iterator mode
    const_iterator begin() const {return const_iterator(this, mMap.begin(), mUseRangeIters);}
    const_iterator cbegin() const {return const_iterator(this, mMap.begin(), mUseRangeIters);}
    const_iterator end() const {return const_iterator(this, mMap.end(), mUseRangeIters);}
    const_iterator cend() const {return const_iterator(this, mMap.end(), mUseRangeIters);}

    // Returns iterators that iterator all values in the span, including non-existent entries that take on the default value
    // The key value is accessible using <iterator>.index()
    // Use <iterator>.element_exists() or <iterator>.element_does_not_exist() to determine whether the entry actually exists in the array or not
    const_iterator contiguous_begin() const {return const_iterator(this, mMap.begin(), true);}
    const_iterator contiguous_cbegin() const {return const_iterator(this, mMap.begin(), true);}
    const_iterator contiguous_end() const {return const_iterator(this, mMap.end(), true);}
    const_iterator contiguous_cend() const {return const_iterator(this, mMap.end(), true);}

    // Returns iterator that iterator only the values entered in the array, excluding non-existent entries that take on the default value
    // The key value is accessible using <iterator>.index()
    // <iterator>.element_exists() should always return true for valid iterators obtained from these functions
    const_iterator skip_begin() const {return const_iterator(this, mMap.begin(), false);}
    const_iterator skip_cbegin() const {return const_iterator(this, mMap.begin(), false);}
    const_iterator skip_end() const {return const_iterator(this, mMap.end(), false);}
    const_iterator skip_cend() const {return const_iterator(this, mMap.end(), false);}

    // Complexity: best case is O(1) (when decayed to a simple vector), worst case is O(logn) (when every element is in its own bucket)
    // If using contiguous iterators, the desired key iterator will be returned if the key is within the current span, end() otherwise
    // If using skip iterators, the desired key iterator will only be returned if a value is found in the array, end() otherwise
    const_iterator iterator_at(Key idx) const {
        return mUseRangeIters? contiguous_iterator_at(idx): skip_iterator_at(idx);
    }
    const_iterator contiguous_iterator_at(Key idx) const {
        auto upper_bound = mMap.upper_bound(idx);
        if (upper_bound == mMap.begin())
            return end();

        --upper_bound;

        if (is_not_in_vector(upper_bound, idx) && upper_bound == --mMap.end())
            return end();

        return const_iterator(this, upper_bound, idx, true);
    }
    const_iterator skip_iterator_at(Key idx) const {
        auto upper_bound = mMap.upper_bound(idx);
        if (upper_bound == mMap.begin())
            return end();

        --upper_bound;

        if (is_not_in_vector(upper_bound, idx))
            return end();

        return const_iterator(this, upper_bound, idx, false);
    }

    // Complexity: O(n)
    void clear() {mMap.clear();}

    // Complexity: best case is O(1) (when decayed to a simple vector), worst case is O(logn) (when every element is in its own bucket)
    const T &operator[](Key idx) const {return at(idx);}
    const T &at(Key idx) const {
        // Find map entry past the desired index
        auto upper_bound = mMap.upper_bound(idx);
        if (upper_bound == mMap.begin())
            return mGenerator(idx);

        // Find map entry the desired index should be in
        --upper_bound;

        // Check if desired index is in bounds
        if (is_not_in_vector(upper_bound, idx))
            return mGenerator(idx);

        // Return element
        return upper_bound->second[idx - upper_bound->first];
    }

    // Complexity: best case is O(1) (when decayed to a simple vector), worst case is O(logn) (when every element is in its own bucket)
    bool contains(Key idx) const {return !skip_iterator_at(idx).atEnd();}

    // Complexity: best case is O(1) (when entire array is decayed to a vector and the element already exists),
    //             average case is O(logn) (when every element is in its own bucket),
    //             worst case is O(n) (when entire array is decayed to a vector and the element doesn't exist)
    T &operator[](Key idx) {return write(idx, mGenerator(idx));}
    T &write(Key idx, const T &item) {
        // Find map entry past the desired index
        auto upper_bound = mMap.upper_bound(idx);
        auto lower_bound = upper_bound;

        // Determine whether the new index requires a vector merge
        if (upper_bound != mMap.begin())
            --lower_bound;

        // If lower_bound == upper_bound, the new item should be the first map entry
        // Insert at start of vector or create a prior map entry
        if (lower_bound == upper_bound) {
            auto it = mMap.insert(upper_bound, std::make_pair(idx, std::vector<T>({item})));
            ++mElements;
            compact(it, upper_bound); // Compact to merge the two vectors together if necessary
            return it->second[0];
        }
        // If lower_bound != upper_bound, the new item will be inserted in a map entry or between map entries
        else {
            // If not in vector, insert as standalone map entry or append to previous
            if (is_not_in_vector(lower_bound, idx)) {
                // Append to end of previous vector?
                if (should_append_to_vector(lower_bound, idx)) {
                    lower_bound->second.push_back(item);
                    ++mElements;
                    compact(lower_bound, upper_bound);
                    return lower_bound->second[idx - vector_begin_idx(lower_bound)];
                }
                // Insert as standalone map entry
                else {
                    auto it = mMap.insert(upper_bound, std::make_pair(idx, std::vector<T>({item})));
                    ++mElements;
                    compact(it, upper_bound);
                    return it->second[0];
                }
            }
            // Otherwise, the desired element is in the vector, so assign and return it
            else
                return lower_bound->second[idx - vector_begin_idx(lower_bound)] = item;
        }
    }

    bool exists(Key key) const {return iterator_at(key).element_exists();}

    // Complexity: best case is O(1) (when decayed to a vector and the specified element is at the end)
    //             average case is O(logn) (when every element is in a bucket by itself)
    //             worst case is O(n) (when decayed to a vector and the specified element is at the beginning)
    void erase(Key key) {
        // Look for element beyond specified element
        auto upper_bound = mMap.upper_bound(key);

        // If the first element, key is not set in the array
        if (upper_bound == mMap.begin())
            return;

        // Get the element that key would be in, if it existed
        --upper_bound;

        // Check if key is in the entry vector
        if (is_not_in_vector(upper_bound, key))
            return;

        // Decrement number of contained elements
        --mElements;

        // If at the end of the vector, just remove it from the vector
        // Falls through to the next case if only one element is in the vector
        if (key == vector_last_used_end_idx(upper_bound) && upper_bound->second.size() > 1) {
            upper_bound->second.pop_back();
        }
        // If at the beginning of the vector, a map rebuild is necessary to fix the start key
        else if (key == vector_begin_idx(upper_bound)) {
            auto it = upper_bound;

            // Only move map entry if there will be something in the vector
            if (upper_bound->second.size() > 1)
            {
                // Move existing vector to new vector
                std::vector<T> vec(std::move(upper_bound->second));
                // Erase beginning of new vector (getting rid of the specified key)
                vec.erase(vec.begin());
                // Insert new vector as a new map entry
                mMap.insert(++upper_bound, std::make_pair(key + 1, std::move(vec)));
            }
            // Erase the existing (now invalid) map entry
            mMap.erase(it);
        }
        // Otherwise, split the vector into two map entries, removing the specified key
        else {
            size_t pivot = key - vector_begin_idx(upper_bound);

            std::vector<T> vec;
            // Copy end of existing vector into new vector
            std::move(upper_bound->second.begin() + pivot + 1, upper_bound->second.end(), std::back_inserter(vec));
            // Erase end of existing vector
            upper_bound->second.erase(upper_bound->second.begin() + pivot, upper_bound->second.end());
            // Insert new vector as a new map entry
            mMap.insert(++upper_bound, std::make_pair(key + 1, std::move(vec)));
        }
    }

    // Complexity: O(1), constant-time
    const T &default_value() const {
        return mDefaultValue;
    }

    // Complexity: O(1), constant-time
    T &set_default_value(const T &value) {
        return mDefaultValue = value;
    }

    // Complexity: O(1), constant-time
    bool empty() const {return mMap.empty();}

    // Complexity: O(1), constant-time
    Key span_begin() const {
        if (empty())
            return Key{};

        return vector_begin_idx(mMap.begin());
    }

    // Complexity: O(1), constant-time
    // The key returned is the last actually used key, not one past the end
    Key span_end() const {
        if (empty())
            return Key{};

        return vector_last_used_end_idx(--mMap.end());
    }

    // Complexity: O(1), constant-time
    // Returns the number of elements this sparse array encompasses (even if they contain the default value)
    // If empty() returns false and this function returns 0, overflow has occured and the span is actually
    // one more than the maximum representable integer value. Keep this in mind when using span() on large key types.
    DiffType span() const {
        if (empty())
            return 0;

        Key first = vector_begin_idx(mMap.begin());
        Key last = vector_last_used_end_idx(--mMap.end());

        return static_cast<DiffType>(last) - static_cast<DiffType>(first) + 1;
    }

    // Complexity: O(1), constant-time
    // Returns true if the given index falls within the currently contained key span
    bool span_contains(Key idx) const {
        if (empty())
            return false;

        return !((idx < vector_begin_idx(mMap.begin())) || vector_last_used_end_idx(--mMap.end()) < idx);
    }

    // Complexity: O(1), constant-time
    // Returns the number of elements in this array
    // Note that though the span can cover the entire range of representable integers (or even more, depending on the key type),
    // the number of storable elements is limited to the range of size_t, no more
    size_t elements() const {return mElements;}
    size_t size() const {return mElements;}

    // Complexity: O(1), constant-time
    // Returns the number of runs (contiguous sections) in this array
    size_t runs() const {return mMap.size();}

    // Complexity: best-cast O(1), worst-case O(n)
    // Returns true if the two sparse arrays are equal, including the default element value
    template<typename K, typename V>
    friend bool operator==(const sparse_array<K, V> &lhs, const sparse_array<K, V> &rhs);

private:
    Key vector_begin_idx(map_const_iterator it) const {
        return it->first;
    }

    Key vector_end_idx(map_const_iterator it) const {
        return it->first + it->second.size();
    }

    Key vector_last_used_end_idx(map_const_iterator it) const {
        return it->first + (it->second.size() - 1);
    }

    bool is_not_in_vector(map_const_iterator it, Key key) const {
        return (key < vector_begin_idx(it) || vector_last_used_end_idx(it) < key);
    }

    bool should_append_to_vector(map_const_iterator it, Key key) const {
        return (vector_last_used_end_idx(it) < key && vector_end_idx(it) == key);
    }

    void compact(map_iterator first,
                 map_iterator second) {
        if (first == second || second == mMap.end() || !(vector_end_idx(first) == vector_begin_idx(second)))
            return;

        first->second.insert(first->second.end(), second->second.begin(), second->second.end());
        mMap.erase(second);
    }
};

namespace impl {
    template<typename Key, typename T, typename DiffType>
    const T &GeneratorHelper<Key, T, DiffType, void>::operator()(Key) {
        return parent->default_value();
    }
}

template<typename Key, typename T>
bool operator==(const sparse_array<Key, T> &lhs, const sparse_array<Key, T> &rhs)
{
    return lhs.mDefaultValue == rhs.mDefaultValue && lhs.mMap == rhs.mMap;
}

template<typename Key, typename T>
bool operator!=(const sparse_array<Key, T> &lhs, const sparse_array<Key, T> &rhs)
{
    return !(lhs == rhs);
}

template<typename Key, typename T, typename DiffType = uintmax_t, typename Generator = void>
sparse_array<Key, T, DiffType, Generator> make_sparse_array(const T &default_value, bool iteratorsEncompassAllValuesInSpan = true) {
    return sparse_array<Key, T, DiffType, Generator>(default_value, iteratorsEncompassAllValuesInSpan);
}

template<typename Key, typename T, typename DiffType = uintmax_t, typename Generator = void>
sparse_array<Key, T, DiffType, Generator> make_sparse_array(const Generator &generator, bool iteratorsEncompassAllValuesInSpan = true) {
    return sparse_array<Key, T, DiffType, Generator>(generator, iteratorsEncompassAllValuesInSpan);
}
