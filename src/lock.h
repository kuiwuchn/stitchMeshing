#pragma once

#include <condition_variable>

class ordered_lock {
public:
    ordered_lock() : next_ticket(0), counter(0) {}
    void lock() {
        std::unique_lock<std::mutex> acquire(cvar_lock);
        unsigned int ticket = next_ticket++;
        while (ticket != counter)
            cvar.wait(acquire);
    }
    void unlock() {
        std::unique_lock<std::mutex> acquire(cvar_lock);
        counter++;
        cvar.notify_all();
    }
protected:
    std::condition_variable  cvar;
    std::mutex               cvar_lock;
    unsigned int             next_ticket, counter;
};

