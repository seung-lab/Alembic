using HttpServer
using WebSockets
using MbedTLS

rel(p::String) = joinpath(dirname(@__FILE__), p)

wsh = WebSocketHandler() do req,client
        while true
            msg = read(client)
            write(client, msg)
        end
      end

cert = MbedTLS.crt_parse_file(rel("./certificate.crt"))
key = MbedTLS.parse_keyfile(rel("./privateKey.key"))

server = Server(wsh)
run(server, port=9999, ssl=(cert, key))