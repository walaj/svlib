/* pthread-lite - Template classes for C++ concurrency
 * Copyright 2016 Jeremiah Wala
 * Modified and expanded by Jeremiah Wala (jwala@broadinstitute.org)
 * Released under the MIT license
 * 
 * Modified from wqueue.h by Vic Hargrave 
 * https://github.com/vichargrave/wqueue
 */

#ifndef PTHREAD_LITE_H__
#define PTHREAD_LITE_H__

#include <list>
#include <pthread.h>

/** Class to hold a queue of work items to
 * be processed by individual threads
 */
template <typename T> class WorkQueue { 

  public:
  WorkQueue() {
    pthread_mutex_init(&m_mutex, NULL);
    pthread_cond_init(&m_condv, NULL);
  }

  ~WorkQueue() {
    pthread_mutex_destroy(&m_mutex);
    pthread_cond_destroy(&m_condv);
  }

  void add(T item) {
    pthread_mutex_lock(&m_mutex);
    m_queue.push_back(item);
    pthread_cond_signal(&m_condv);
    pthread_mutex_unlock(&m_mutex);
  }

  T remove() {
    pthread_mutex_lock(&m_mutex);
    while (m_queue.size() == 0) {
      pthread_cond_wait(&m_condv, &m_mutex);
    }
    T item = m_queue.front();
    m_queue.pop_front();
    pthread_mutex_unlock(&m_mutex);
    return item;
  }

  int size() {
    pthread_mutex_lock(&m_mutex);
    int size = m_queue.size();
    pthread_mutex_unlock(&m_mutex);
    return size;
  }

  std::list<T>   m_queue;
  pthread_mutex_t m_mutex;
  pthread_cond_t  m_condv;

};

class WorkThread {

  public:
  WorkThread() : m_tid(0), m_running(0), m_detached(0) {}
  
  ~WorkThread()
  {
    if (m_running == 1 && m_detached == 0) {
      pthread_detach(m_tid);
    }
    if (m_running == 1) {
      pthread_cancel(m_tid);
    }
  }

  static void* runThread(void* arg) {
    return ((WorkThread*)arg)->run();
  }
 
  int start()
  {
    int result = pthread_create(&m_tid, NULL, runThread, this);
    if (result == 0) {
      m_running = 1;
    }
    return result;
  }

  int join()
  {
    int result = -1;
    if (m_running == 1) {
      result = pthread_join(m_tid, NULL);
      if (result == 0) {
	m_detached = 1;
      }
    }
    return result;
  }
 
  int detach()
  {
    int result = -1;
    if (m_running == 1 && m_detached == 0) {
      result = pthread_detach(m_tid);
      if (result == 0) {
	m_detached = 1;
      }
    }
    return result;
  }

  pthread_t self() {
    return m_tid;
  }

  virtual void* run() = 0;
 
private:
  pthread_t  m_tid;
  int        m_running;
  int        m_detached;
};


/** Class to hold both a reference to a queue of
 * work items, and a set of data specific to a particular
 * thread (eg separate pointer to a location in a file
 *  to avoid thread collisions when randomly accessing files)
 */
template <typename W, typename T> // W work item, T thread item
class ConsumerThread : public WorkThread {
public:

 ConsumerThread(WorkQueue<W*>& queue, T* thread_data) : m_queue(queue), m_thread_data(thread_data) {}
  
  void* run() {
    // Remove 1 item at a time and process it. Blocks if no items are 
    // available to process.
    for (int i = 0;; i++) {
      W* item = (W*)m_queue.remove();
      if (!item)
	return NULL;
      item->run(m_thread_data); // pass thread data to worker unit always
      delete item;
      if (m_queue.size() == 0)
        return NULL;
    }
    return NULL;
  }

  // return thread data. Useful at end when want to process
  // results stored per each thread
  T* GetThreadData() const { return m_thread_data; }

 private: 
  WorkQueue<W*>& m_queue; // a queue of worker units
  T* m_thread_data; // a set of data private to this thread (often not uses)
  
};


#endif
