FROM nixos/nix:latest as builder

RUN mkdir -p /etc/nix && \
    echo "experimental-features = nix-command flakes" >> /etc/nix/nix.conf

WORKDIR /app

COPY . .

RUN nix build .

RUN mkdir nix-store-closure

RUN cp -R $(nix-store -qR result/) nix-store-closure

FROM scratch

WORKDIR /app

# Copy /nix/store
COPY --from=builder /app/nix-store-closure /nix/store
COPY --from=builder /app/result /app
CMD ["/app/bin/highs"]
