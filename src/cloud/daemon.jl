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
        message = Queue.pop_message(daemon.queue)

        task = DaemonTask.parse(message)

        print("Task is $(task.taskId)")

        #=DaemonTask.execute(task)=#

        sleep(poll_frequency_seconds)
    end
end

end # end module Daemon
