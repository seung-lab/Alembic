module Daemon

import Julimaps.Cloud.Queue
import Julimaps.Cloud.Bucket
import Julimaps.Cloud.DaemonTask

export DaemonService
export run

type DaemonService
    queue::Queue.QueueService
    bucket::Bucket.BucketService
    poll_frequency_seconds::Int64
end

function run(daemon::DaemonService)
    while true
        try
            message = Queue.pop_message(daemon.queue)

            task = DaemonTask.parse(message)

            print("Task is $(task.taskId)")

            #=DaemonTask.execute(task)=#
        catch e
            showerror(STDERR, e, catch_backtrace(); backtrace = true)
            print(STDERR, "\n") #looks like showerror doesn't include a newline
        end

        sleep(daemon.poll_frequency_seconds)
    end
end

end # end module Daemon
