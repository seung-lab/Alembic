module Daemon
import AWS
import QueueService
import BucketService
import DaemonTask

type Daemon
    queue::QueueService.QueueService
    bucket::BucketService.BucketService
    poll_frequency_seconds::Int64
end

function run(daemon::Daemon)
    while true
        message = QueueService.pop_message(daemon.queue)
        task = DaemonTask.parse(message)
        print("Task is $(task.taskId)")

        DaemonTask.execute(task)
        sleep(poll_frequency_seconds)
    end
end

end # end module Daemon
