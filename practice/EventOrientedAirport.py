# Event oriented simulation of airport
# In class example


import Queue


# Globals
fel = None
in_the_air = 0
on_the_ground = 0
runway_free = True
r = 10
g = 60


def schedule_event(event):
    global fel
    fel.put(event)


def arrival_event(now):
    global in_the_air, runway_free
    in_the_air += 1
    if runway_free:
        runway_free = False
        schedule_event((now + r, 1))
    display_state(now)


def landed_event(now):
    global on_the_ground, in_the_air, runway_free
    in_the_air -= 1
    on_the_ground += 1
    schedule_event((now + g, 2))
    runway_free = True
    display_state(now)


def departure_event(now):
    global on_the_ground
    on_the_ground -= 1
    display_state(now)


def display_state(time_stamp):
    global on_the_ground, in_the_air, runway_free
    print(time_stamp, in_the_air, on_the_ground, runway_free)


if __name__ == "__main__":

    # Initial List of events
    fel = Queue.PriorityQueue(maxsize=100)
    fel.put((0, 0))   # Landed event at 0 min
    fel.put((20, 0))  # Landed event at 20 min

    # Events are represented as tuples (time, type), where time is the event's
    # scheduled time and type is the event type (0 = arrive, 1 = landed, and
    # 2 = departure). The time stamp is the first element in the tuple so that
    # the earliest event will be popped from the priority queue.

    # Event processing loop
    while not fel.empty():
        e = fel.get()
        n = e[0]
        if e[1] == 0:
            arrival_event(n)
        elif e[1] == 1:
            landed_event(n)
        elif e[1] == 2:
            departure_event(n)
        else:
            print("Invalid event")
