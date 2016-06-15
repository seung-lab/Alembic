module Daemon

import Julimaps.Cloud.Queue
import Julimaps.Cloud.Bucket

export DaemonService
export run

type DaemonService
    queue::Queue.AWSQueueService
    bucket::Bucket.AWSBucketService
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
